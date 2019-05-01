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
from microhapulator.genotype import ObservedGenotype


def dist(gt1, gt2):
    allloci = set(gt1.alleles()).union(gt2.alleles())
    hammdist = 0
    for locus in allloci:
        allele1, allele2 = None, None
        if locus in gt1.data:
            allele1 = gt1.data[locus]['genotype']
        if locus in gt2.data:
            allele2 = gt2.data[locus]['genotype']
        if allele1 != allele2:
            hammdist += 1
    return hammdist


def main(args):
    gt1 = ObservedGenotype(filename=args.gt1)
    gt2 = ObservedGenotype(filename=args.gt2)
    d = dist(gt1, gt2)
    data = {
        'hamming_distance': d,
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
