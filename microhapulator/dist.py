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
from microhapulator.genotype import Genotype


def dist(gt1, gt2):
    allloci = set(gt1.loci()).union(gt2.loci())
    hammdist = 0
    for locus in allloci:
        allele1 = gt1.alleles(locus)
        allele2 = gt2.alleles(locus)
        if allele1 != allele2:
            hammdist += 1
    return hammdist


def main(args):
    d = dist(Genotype(fromfile=args.genotype1), Genotype(fromfile=args.genotype2))
    data = {
        'hamming_distance': d,
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
