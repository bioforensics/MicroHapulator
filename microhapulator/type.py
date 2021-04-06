#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
import microhapulator
import os.path
import pysam
import re


def parse_variant_offsets_from_fasta_headers(fasta):
    offsets = dict()
    for line in fasta:
        if not line.startswith('>'):
            continue
        markerid = line[1:].split()[0]
        if ' variants=' not in line:
            message = 'variant offsets not annotated for target amplicon: ' + line
            raise ValueError(message)
        offsetstr = re.search(r'variants=(\S+)', line).group(1)
        varloc = [int(x) for x in offsetstr.split(',')]
        offsets[markerid] = varloc
    return offsets


def check_index(bamfile):
    index1 = bamfile + '.bai'
    index2 = re.sub(r'.bam$', '.bai', bamfile)
    for testfile in (index1, index2):
        if os.path.isfile(testfile):
            break
    else:
        pysam.index(bamfile)


def tally_haplotypes(bamfile, refrfasta):
    totaldiscarded = 0
    with microhapulator.open(refrfasta, 'r') as fh:
        offsets = parse_variant_offsets_from_fasta_headers(fh)
    check_index(bamfile)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    for locusid in sorted(offsets):
        discarded = 0
        haplotypes = defaultdict(int)
        ht = defaultdict(dict)
        varloc = set(offsets[locusid])
        cov_pos = list()
        for column in bam.pileup(locusid):
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
            htstr = ','.join(htlist)
            haplotypes[htstr] += 1
        yield locusid, cov_pos, haplotypes, discarded
        totaldiscarded += discarded
    microhapulator.plog(
        '[MicroHapulator::type] discarded', totaldiscarded,
        'reads with gaps or missing data at positions of interest'
    )


def type(bamfile, refrfasta, threshold=None):
    genotyper = tally_haplotypes(bamfile, refrfasta)
    profile = microhapulator.profile.ObservedProfile()
    for locusid, cov_by_pos, htcounts, ndiscarded in genotyper:
        profile.record_coverage(locusid, cov_by_pos, ndiscarded=ndiscarded)
        for allele, count in htcounts.items():
            profile.record_allele(locusid, allele, count)
        if threshold is None:
            profile.data['markers'][locusid]['genotype'] = list()
    if threshold is not None:
        profile.infer(threshold=threshold)
    return profile


def main(args):
    profile = type(args.bam, args.refr, threshold=args.threshold)
    profile.dump(args.out)
