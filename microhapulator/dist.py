# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

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
        "hamming_distance": d,
    }
    with microhapulator.open(args.out, "w") as fh:
        json.dump(data, fh, indent=4)
