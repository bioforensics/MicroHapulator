#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

# Core library imports
from collections import namedtuple
from os import fsync
from shutil import rmtree
from string import ascii_letters, digits
from subprocess import check_call
import sys
from tempfile import mkdtemp

# Third-party library imports
from happer.mutate import mutate
from numpy.random import choice, randint

# Internal imports
import microhapulator
from microhapulator.profile import Profile


SimulatedRead = namedtuple("SimulatedRead", ["identifier", "sequence", "quality"])


def calc_n_reads_from_proportions(n, totalreads, prop):
    if prop is None:
        prop = [1.0 / n for _ in range(n)]
    else:
        if len(prop) != n:
            raise ValueError("mismatch between contributor number and proportions")
    normprop = [x / sum(prop) for x in prop]
    return [int(totalreads * x) for x in normprop]


def new_signature():
    return "".join([choice(list(ascii_letters + digits)) for _ in range(7)])


def sequencing(profile, seed=None, threads=1, numreads=500000, readsignature=None, readindex=0):
    tempdir = mkdtemp()
    try:
        haplofile = tempdir + "/haplo.fasta"
        with microhapulator.open(haplofile, "w") as fh:
            for defline, sequence in profile.haploseqs:
                print(">", defline, "\n", sequence, sep="", file=fh)
        isscmd = [
            "iss",
            "generate",
            "--n_reads",
            str(numreads),
            "--draft",
            haplofile,
            "--model",
            "MiSeq",
            "--output",
            tempdir + "/seq",
            "--quiet",
        ]
        if seed:
            isscmd.extend(["--seed", str(seed)])
        if threads:
            isscmd.extend(["--cpus", str(threads)])
        try:
            microhapulator.logstream.flush()
            fsync(microhapulator.logstream.fileno())
        except (AttributeError, OSError):  # pragma: no cover
            pass
        check_call(isscmd)
        f1, f2 = tempdir + "/seq_R1.fastq", tempdir + "/seq_R2.fastq"
        with open(f1, "r") as infh1, open(f2, "r") as infh2:
            if readsignature is None:
                readsignature = new_signature()
            linebuffer = list()
            for line_r1, line_r2 in zip(infh1, infh2):
                if line_r1.startswith("@mh"):
                    readindex += 1
                    prefix = f"@{readsignature}:0:0:0:0:0:{readindex} 1:N:0:0 mh"
                    line_r1 = line_r1.replace("@mh", prefix, 1)
                    prefix = f"@{readsignature}:0:0:0:0:0:{readindex} 2:N:0:0 mh"
                    line_r2 = line_r2.replace("@mh", prefix, 1)
                    line_r1 = line_r1[:-3] + "\n"
                    line_r2 = line_r2[:-3] + "\n"
                linebuffer.append((line_r1, line_r2))
                if len(linebuffer) == 4:
                    r1 = SimulatedRead(
                        identifier=linebuffer[0][0],
                        sequence=linebuffer[1][0],
                        quality=linebuffer[3][0],
                    )
                    r2 = SimulatedRead(
                        identifier=linebuffer[0][1],
                        sequence=linebuffer[1][1],
                        quality=linebuffer[3][1],
                    )
                    yield readindex, r1, r2
                    linebuffer = list()
    finally:
        rmtree(tempdir)


def seq(profiles, seeds=None, threads=1, totalreads=500000, proportions=None, sig=None):
    n = len(profiles)
    if seeds is None:
        seeds = [randint(1, 2 ** 32 - 1) for _ in range(n)]
    if len(seeds) != n:
        raise ValueError("number of profiles must match number of seeds")
    numreads = calc_n_reads_from_proportions(n, totalreads, proportions)
    if 0 in numreads:
        raise ValueError("specified proportions result in 0 reads for 1 or more individuals")
    readsignature = sig if sig else new_signature()
    reads_sequenced = 0
    for profile, seed, nreads in zip(profiles, seeds, numreads):
        message = "Individual seed={seed} numreads={n}".format(seed=seed, n=nreads)
        microhapulator.plog("[MicroHapulator::seq]", message)
        sequencer = sequencing(
            profile,
            seed=seed,
            threads=threads,
            numreads=nreads,
            readsignature=readsignature,
            readindex=reads_sequenced,
        )
        for data in sequencer:
            yield data
        reads_sequenced = data[0]


def resolve_profiles(gtfiles):
    profiles = list()
    for gtfile in gtfiles:
        profile = Profile(fromfile=gtfile)
        for p in profile.unmix():
            profiles.append(p)
    return profiles


def main(args):
    if len(args.out) not in (1, 2):
        raise ValueError(f"expected 1 or 2 output files, found {len(args.out)}")
    fh1 = fh2 = args.out[0]
    if len(args.out) == 2:
        fh2 = args.out[1]
    profiles = resolve_profiles(args.profiles)
    sequencer = seq(
        profiles,
        seeds=args.seeds,
        threads=args.threads,
        totalreads=args.num_reads,
        proportions=args.proportions,
        sig=args.signature,
    )
    for n, read1, read2 in sequencer:
        print(read1.identifier, read1.sequence, "+\n", read1.quality, sep="", end="", file=fh1)
        print(read2.identifier, read2.sequence, "+\n", read2.quality, sep="", end="", file=fh2)
    for fh in args.out:
        if fh != sys.stdout:
            fh.close()
