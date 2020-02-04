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


def diff(prof1, prof2):
    allmarkers = set(prof1.markers()).union(prof2.markers())
    for marker in sorted(allmarkers):
        allele1 = prof1.alleles(marker)
        allele2 = prof2.alleles(marker)
        diff1 = allele1 - allele2
        diff2 = allele2 - allele1
        if len(diff1) > 0 or len(diff2) > 0:
            yield marker, diff1, diff2


def main(args):
    differ = diff(Profile(fromfile=args.profile1), Profile(fromfile=args.profile2))
    with microhapulator.open(args.out, 'w') as fh:
        for marker, diff1, diff2 in differ:
            print(marker, file=fh)
            if len(diff1) > 0:
                for haplotype in sorted(diff1):
                    print('>>>', haplotype, file=fh)
            if len(diff2) > 0:
                for haplotype in sorted(diff2):
                    print('<<<', haplotype, file=fh)
