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


def dist(p1, p2):
    hammdist = 0
    for marker in set(p1.markers()).union(p2.markers()):
        allele1 = p1.alleles(marker)
        allele2 = p2.alleles(marker)
        if allele1 != allele2:
            hammdist += 1
    return hammdist


def main(args):
    d = dist(Profile(fromfile=args.profile1), Profile(fromfile=args.profile2))
    data = {
        'hamming_distance': d,
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
