#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
import microhapulator
from microhapulator.genotype import ObservedGenotype, SimulatedGenotype


def dist(gt1, gt2):
    allloci = set(gt1.loci()).union(gt2.loci())
    hammdist = 0
    for locus in allloci:
        allele1, allele2 = None, None
        if locus in gt1.data:
            allele1 = gt1.alleles(locus)
        if locus in gt2.data:
            allele2 = gt2.alleles(locus)
        if allele1 != allele2:
            hammdist += 1
    return hammdist


def main(args):
    if len(args.sim) + len(args.obs) != 2:
        raise ValueError('the number of BED and/or JSON files must be 2!')
    genotypes = list()
    for bedfile in args.sim:
        with microhapulator.open(bedfile, 'r') as fh:
            gt = SimulatedGenotype(frombed=fh)
        genotypes.append(gt)
    for jsonfile in args.obs:
        gt = ObservedGenotype(filename=jsonfile)
        genotypes.append(gt)
    d = dist(*genotypes)
    data = {
        'hamming_distance': d,
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
