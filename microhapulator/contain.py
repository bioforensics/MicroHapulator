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
from microhapulator.profile import Profile


def contain(gt1, gt2):
    '''Compute the proportion of alleles from gt2 present in gt1.'''
    total = 0
    contained = 0
    for locus in gt2.loci():
        allele1 = gt1.alleles(locus)
        allele2 = gt2.alleles(locus)
        total += len(allele2)
        contained += len(allele2 & allele1)
    return contained, total


def main(args):
    contained, total = contain(
        Genotype(fromfile=args.genotype1),
        Genotype(fromfile=args.genotype2)
    )
    data = {
        'containment': round(contained / total, 4),
        'contained_alleles': contained,
        'total_alleles': total
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
