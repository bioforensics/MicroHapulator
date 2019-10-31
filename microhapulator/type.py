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


class MissingBAMIndexError(ValueError):
    pass


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
        message = 'Please index "{bam:s}" with "samtools index"'.format(bam=bamfile)
        raise MissingBAMIndexError(message)


def observe_genotypes(bamfile, refrfasta):
    totaldiscarded = 0
    with microhapulator.open(refrfasta, 'r') as fh:
        offsets = parse_variant_offsets_from_fasta_headers(fh)
    check_index(bamfile)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    for locusid in sorted(offsets):
        discarded = 0
        genotypes = defaultdict(int)
        gt = defaultdict(dict)
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
                gt[record.alignment.query_name][column.pos] = aligned_base
        for readname, gtdict in gt.items():
            gtlist = [gtdict[pos] for pos in sorted(gtdict)]
            if len(gtlist) < len(varloc):
                discarded += 1
                continue
            gtstr = ','.join(gtlist)
            genotypes[gtstr] += 1
        yield locusid, cov_pos, genotypes, discarded
        totaldiscarded += discarded
    microhapulator.plog(
        '[MicroHapulator::type] discarded', totaldiscarded,
        'reads with gaps or missing data at positions of interest'
    )


def type(bamfile, refrfasta, threshold=10):
    genotyper = observe_genotypes(bamfile, refrfasta)
    gt = microhapulator.profile.ObservedProfile()
    for locusid, cov_by_pos, gtcounts, ndiscarded in genotyper:
        gt.record_coverage(locusid, cov_by_pos, ndiscarded=ndiscarded)
        for allele, count in gtcounts.items():
            gt.record_allele(locusid, allele, count)
    gt.infer(threshold=threshold)
    return gt


def main(args):
    gt = type(args.bam, args.refr, threshold=args.threshold)
    gt.dump(args.out)
