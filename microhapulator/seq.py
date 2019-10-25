#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

# Core library imports
from os import fsync
from shutil import rmtree
from string import ascii_letters, digits
from subprocess import check_call
from tempfile import mkdtemp

# Third-party library imports
from happer.mutate import mutate
from numpy.random import choice, randint

# Internal imports
import microhapulator
from microhapulator.profile import Profile


def calc_n_reads_from_proportions(n, totalreads, prop):
    if prop is None:
        prop = [1.0 / n for _ in range(n)]
    else:
        if len(prop) != n:
            raise ValueError('mismatch between contributor number and proportions')
    normprop = [x / sum(prop) for x in prop]
    return [int(totalreads * x) for x in normprop]


def new_signature():
    return ''.join([choice(list(ascii_letters + digits)) for _ in range(7)])


def sequencing(genotype, seed=None, threads=1, numreads=500000,
               readsignature=None, readindex=0, debug=False):
    tempdir = mkdtemp()
    try:
        haplofile = tempdir + '/haplo.fasta'
        with microhapulator.open(haplofile, 'w') as fh:
            for defline, sequence in genotype.haploseqs:
                print('>', defline, '\n', sequence, sep='', file=fh)
        isscmd = [
            'iss', 'generate', '--n_reads', str(numreads * 2), '--draft', haplofile,
            '--model', 'MiSeq', '--output', tempdir + '/seq', '--quiet'
        ]
        if seed:
            isscmd.extend(['--seed', str(seed)])
        if threads:
            isscmd.extend(['--cpus', str(threads)])
        try:
            microhapulator.logstream.flush()
            fsync(microhapulator.logstream.fileno())
        except (AttributeError, OSError):  # pragma: no cover
            pass
        if debug:
            check_call(isscmd, stderr=microhapulator.logstream)
        else:
            check_call(isscmd)
        with open(tempdir + '/seq_R1.fastq', 'r') as infh:
            if readsignature is None:
                readsignature = new_signature()
            linebuffer = list()
            for line in infh:
                if line.startswith('@MHDBL'):
                    readindex += 1
                    prefix = '@{sig:s}_read{n:d} MHDBL'.format(sig=readsignature, n=readindex)
                    line = line.replace('@MHDBL', prefix, 1)
                linebuffer.append(line)
                if len(linebuffer) == 4:
                    yield readindex, linebuffer[0], linebuffer[1], linebuffer[3]
                    linebuffer = list()
    finally:
        rmtree(tempdir)


def seq(genotypes, seeds=None, threads=1, totalreads=500000, proportions=None,
        sig=None, debug=False):
    n = len(genotypes)
    if seeds is None:
        seeds = [randint(1, 2**32 - 1) for _ in range(n)]
    if len(seeds) != n:
        raise ValueError('number of genotypes must match number of seeds')
    numreads = calc_n_reads_from_proportions(n, totalreads, proportions)
    if 0 in numreads:
        raise ValueError('specified proportions result in 0 reads for 1 or more individuals')
    readsignature = sig if sig else new_signature()
    reads_sequenced = 0
    for genotype, seed, nreads in zip(genotypes, seeds, numreads):
        message = 'Individual seed={seed} numreads={n}'.format(seed=seed, n=nreads)
        microhapulator.plog('[MicroHapulator::seq]', message)
        sequencer = sequencing(
            genotype, seed=seed, threads=threads, numreads=nreads,
            readsignature=readsignature, readindex=reads_sequenced,
            debug=debug,
        )
        for data in sequencer:
            yield data
        reads_sequenced = data[0]


def resolve_genotypes(gtfiles):
    genotypes = list()
    for gtfile in gtfiles:
        genotype = Genotype(fromfile=gtfile)
        for gt in genotype.unmix():
            genotypes.append(gt)
    return genotypes


def main(args):
    genotypes = resolve_genotypes(args.genotypes)
    sequencer = seq(
        genotypes, seeds=args.seeds, threads=args.threads, totalreads=args.num_reads,
        proportions=args.proportions, debug=args.debug, sig=args.signature
    )
    with microhapulator.open(args.out, 'w') as fh:
        for n, defline, sequence, qualities in sequencer:
            print(defline, sequence, '+\n', qualities, sep='', end='', file=fh)
