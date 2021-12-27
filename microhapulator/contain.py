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


def contain(p1, p2):
    """Compute the proportion of alleles from p2 present in p1."""
    total = 0
    contained = 0
    for marker in p2.markers():
        allele1 = p1.alleles(marker)
        allele2 = p2.alleles(marker)
        total += len(allele2)
        contained += len(allele2 & allele1)
    return contained, total


def main(args):
    contained, total = contain(Profile(fromfile=args.profile1), Profile(fromfile=args.profile2))
    data = {
        "containment": round(contained / total, 4),
        "contained_alleles": contained,
        "total_alleles": total,
    }
    with microhapulator.open(args.out, "w") as fh:
        json.dump(data, fh, indent=4)
