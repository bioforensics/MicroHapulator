# -------------------------------------------------------------------------------------------------
# Copyright (c) 2021, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


from collections import namedtuple, defaultdict
from math import ceil
from microhapulator.parsers import load_marker_definitions, cross_check_marker_ids
from microhapulator.profile import SimulatedProfile, TypingResult
import numpy as np
import os
import pandas as pd
import pysam
import re
from shutil import rmtree
from string import ascii_letters, digits
from subprocess import check_call, run
import sys
from tempfile import mkdtemp, TemporaryDirectory


SimulatedRead = namedtuple("SimulatedRead", ["identifier", "sequence", "quality"])


def count_and_sort(profile, include_discarded=True):
    counts = dict(
        Marker=list(),
        ReadCount=list(),
    )
    for marker, mdata in profile.data["markers"].items():
        readcount = 0
        if include_discarded:
            readcount += mdata["num_discarded_reads"]
        for haplotype, count in mdata["allele_counts"].items():
            readcount += count
        counts["Marker"].append(marker)
        counts["ReadCount"].append(readcount)
    data = pd.DataFrame(counts).sort_values(["ReadCount"], ascending=False).reset_index(drop=True)
    return data


def balance(profile, include_discarded=True):
    data = count_and_sort(profile, include_discarded=include_discarded)
    with TemporaryDirectory() as tempdir:
        tfile = os.path.join(tempdir, "data.tsv")
        data.to_csv(tfile, index=False, header=False)
        run(["termgraph", tfile])
    return data


def contain(p1, p2):
    """Compute the proportion of alleles from p2 present in p1."""
    total = 0
    contained = 0
    for marker in p2.markers():
        allele1 = p1.alleles(marker)
        allele2 = p2.alleles(marker)
        total += len(allele2)
        contained += len(allele2 & allele1)
    return contained, total


def contrib(profile):
    num_alleles_per_marker = [len(profile.alleles(marker)) for marker in profile.markers()]
    max_num_alleles = max(num_alleles_per_marker)
    max_thresh = max_num_alleles - 1 if max_num_alleles % 2 == 0 else max_num_alleles
    max_loci = sum([1 for n in num_alleles_per_marker if n >= max_thresh])
    max_perc = round(max_loci / len(num_alleles_per_marker), 4)
    return ceil(max_num_alleles / 2), max_loci, max_perc


def diff(prof1, prof2):
    allmarkers = set(prof1.markers()).union(prof2.markers())
    for marker in sorted(allmarkers):
        allele1 = prof1.alleles(marker)
        allele2 = prof2.alleles(marker)
        diff1 = allele1 - allele2
        diff2 = allele2 - allele1
        if len(diff1) > 0 or len(diff2) > 0:
            yield marker, diff1, diff2


def dist(p1, p2):
    hammdist = 0
    for marker in set(p1.markers()).union(p2.markers()):
        allele1 = p1.alleles(marker)
        allele2 = p2.alleles(marker)
        if allele1 != allele2:
            hammdist += 1
    return hammdist


def prob(frequencies, prof1, prof2=None, erate=0.001):
    if prof2 is None:
        return prof1.rand_match_prob(frequencies)
    else:
        return prof1.rmp_lr_test(prof2, frequencies, erate=erate)


def calc_n_reads_from_proportions(n, totalreads, prop):
    if prop is None:
        prop = [1.0 / n for _ in range(n)]
    else:
        if len(prop) != n:
            raise ValueError("mismatch between contributor number and proportions")
    normprop = [x / sum(prop) for x in prop]
    return [int(totalreads * x) for x in normprop]


def new_signature():
    return "".join([np.random.choice(list(ascii_letters + digits)) for _ in range(7)])


def sequencing(
    profile,
    markers,
    refrseqs,
    seed=None,
    threads=1,
    numreads=500000,
    readsignature=None,
    readindex=0,
):
    tempdir = mkdtemp()
    try:
        haplofile = tempdir + "/haplo.fasta"
        with open(haplofile, "w") as fh:
            for defline, sequence in profile.haploseqs(markers, refrseqs):
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


def seq(
    profiles,
    markers,
    refrseqs,
    seeds=None,
    threads=1,
    totalreads=500000,
    proportions=None,
    sig=None,
):
    n = len(profiles)
    if seeds is None:
        seeds = [np.random.randint(1, 2 ** 32 - 1) for _ in range(n)]
    if len(seeds) != n:
        raise ValueError("number of profiles must match number of seeds")
    numreads = calc_n_reads_from_proportions(n, totalreads, proportions)
    if 0 in numreads:
        raise ValueError("specified proportions result in 0 reads for 1 or more individuals")
    readsignature = sig if sig else new_signature()
    reads_sequenced = 0
    for profile, seed, nreads in zip(profiles, seeds, numreads):
        message = "Individual seed={seed} numreads={n}".format(seed=seed, n=nreads)
        print("[MicroHapulator::seq]", message, file=sys.stderr)
        sequencer = sequencing(
            profile,
            markers,
            refrseqs,
            seed=seed,
            threads=threads,
            numreads=nreads,
            readsignature=readsignature,
            readindex=reads_sequenced,
        )
        for data in sequencer:
            yield data
        reads_sequenced = data[0]


def sim(frequencies, seed=None):
    """Simulate a diploid genotype from the specified microhaplotype frequencies."""
    profile = SimulatedProfile(ploidy=2)
    if seed is None:
        seed = np.random.randint(2 ** 32 - 1)
    profile.data["metadata"] = {
        "HaploSeed": seed,
    }
    np.random.seed(seed)
    markers = sorted(frequencies.Marker.unique())
    for haploindex in range(2):
        for marker in markers:
            haplofreqs = frequencies[frequencies.Marker == marker]
            haplotypes = list(haplofreqs.Haplotype)
            freqs = list(haplofreqs.Frequency)
            freqs = [x / sum(freqs) for x in freqs]
            sampled_haplotype = np.random.choice(haplotypes, p=freqs)
            profile.add(haploindex, marker, sampled_haplotype)
    message = f"simulated microhaplotype variation at {len(markers)} markers"
    print("[MicroHapulator::sim]", message, file=sys.stderr)
    return profile


def check_index(bamfile):
    index1 = bamfile + ".bai"
    index2 = re.sub(r".bam$", ".bai", bamfile)
    for testfile in (index1, index2):
        if os.path.isfile(testfile):
            break
    else:
        pysam.index(bamfile)


def tally_haplotypes(bam, offsets, minbasequal=10, max_depth=1e6):
    totaldiscarded = 0
    for locusid in sorted(offsets):
        discarded = 0
        haplotypes = defaultdict(int)
        ht = defaultdict(dict)
        varloc = set(offsets[locusid])
        cov_pos = list()
        for column in bam.pileup(locusid, min_base_quality=minbasequal, max_depth=max_depth):
            cov_pos.append(column.n)
            if column.pos not in varloc:
                continue
            for record in column.pileups:
                aligned_base = None
                if record.is_del or record.is_refskip:
                    continue
                aligned_base = record.alignment.query_sequence[record.query_position]
                ht[record.alignment.query_name][column.pos] = aligned_base
        for readname, htdict in ht.items():
            htlist = [htdict[pos] for pos in sorted(htdict)]
            if len(htlist) < len(varloc):
                discarded += 1
                continue
            htstr = ",".join(htlist)
            haplotypes[htstr] += 1
        yield locusid, cov_pos, haplotypes, discarded
        totaldiscarded += discarded
    print(
        "[MicroHapulator::type] discarded",
        totaldiscarded,
        "reads with gaps or missing data at positions of interest",
        file=sys.stderr,
    )


def type(
    bamfile, markertsv, minbasequal=10, ecthreshold=0.25, static=None, dynamic=None, max_depth=1e6
):
    check_index(bamfile)
    bam = pysam.AlignmentFile(bamfile, "rb")
    markers = load_marker_definitions(markertsv)
    offsets = defaultdict(list)
    for n, row in markers.iterrows():
        offsets[row.Marker].append(row.Offset)
    cross_check_marker_ids(bam.references, offsets.keys(), "read alignments", "marker definitions")
    genotyper = tally_haplotypes(bam, offsets, minbasequal=minbasequal, max_depth=max_depth)
    result = TypingResult()
    for locusid, cov_by_pos, htcounts, ndiscarded in genotyper:
        result.record_coverage(locusid, cov_by_pos, ndiscarded=ndiscarded)
        for allele, count in htcounts.items():
            result.record_allele(locusid, allele, count)
    result.infer(ecthreshold=ecthreshold, static=static, dynamic=dynamic)
    return result