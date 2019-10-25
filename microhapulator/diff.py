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


def diff(gt1, gt2):
    allloci = set(gt1.loci()).union(gt2.loci())
    for locus in sorted(allloci):
        allele1 = gt1.alleles(locus)
        allele2 = gt2.alleles(locus)
        diff1 = allele1 - allele2
        diff2 = allele2 - allele1
        if len(diff1) > 0 or len(diff2) > 0:
            yield locus, diff1, diff2


def main(args):
    differ = diff(Genotype(fromfile=args.genotype1), Genotype(fromfile=args.genotype2))
    with microhapulator.open(args.out, 'w') as fh:
        for locus, diff1, diff2 in differ:
            print(locus, file=fh)
            if len(diff1) > 0:
                print('>>>', ':'.join(sorted(diff1)), file=fh)
            if len(diff2) > 0:
                print('<<<', ':'.join(sorted(diff2)), file=fh)
