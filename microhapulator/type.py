#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
import json
import microhapulator
import pysam
import sys


class Genotype(object):
    def __init__(self, filename=None):
        self.data = dict()
        if filename:
            with microhapulator.open(filename, 'r') as fh:
                self.data = json.load(fh)

    def record_coverage(self, locusid, cov_by_pos, ndiscarded=0):
        self.data[locusid] = {
            'mean_coverage': round(sum(cov_by_pos) / len(cov_by_pos), 1),
            'min_coverage': min(cov_by_pos),
            'max_coverage': max(cov_by_pos),
            'num_discarded_reads': ndiscarded,
            'allele_counts': dict(),
        }

    def record_allele(self, locusid, allele, count):
        self.data[locusid]['allele_counts'][allele] = count

    def infer(self):
        for locusid, locusdata in self.data.items():
            allelecounts = locusdata['allele_counts']
            avgcount = sum(allelecounts.values()) / len(allelecounts.values())
            gt = set()
            for allele, count in allelecounts.items():
                if count * 2 < avgcount:
                    continue
                gt.add(allele)
            self.data[locusid]['genotype'] = sorted(gt)

    def dump(self, file=None):
        if file is None:
            return json.dumps(self.data, indent=4, sort_keys=True)
        else:
            return json.dump(self.data, file, indent=4, sort_keys=True)


def parse_variant_offsets_from_fasta_headers(fasta):
    offsets = dict()
    for line in fasta:
        if not line.startswith('>'):
            continue
        locusid, refrloc, varinfo = line[1:].strip().split()
        varoffsets = varinfo.split('=')[1]
        varloc = [int(x) for x in varoffsets.split(':')]
        offsets[locusid] = varloc
    return offsets


def observe_genotypes(bamfile, refrfasta):
    totaldiscarded = 0
    with microhapulator.open(refrfasta, 'r') as fh:
        offsets = parse_variant_offsets_from_fasta_headers(fh)
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


def genotype(bamfile, refrfasta):
    genotyper = observe_genotypes(bamfile, refrfasta)
    gt = Genotype()
    for locusid, cov_by_pos, gtcounts, ndiscarded in genotyper:
        gt.record_coverage(locusid, cov_by_pos, ndiscarded=ndiscarded)
        for allele, count in gtcounts.items():
            gt.record_allele(locusid, allele, count)
    gt.infer()
    return gt


def main(args):
    gt = genotype(args.bam, args.refr)
    with microhapulator.open(args.out, 'w') as fh:
        gt.dump(file=fh)
