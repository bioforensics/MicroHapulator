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
from microhapulator.parsers import open as mhopen
import microhapulator.api as mhapi
from microhapulator.profile import Profile


def subparser(subparsers):
    cli = subparsers.add_parser("contain")
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument("profile1", help="simulated or inferred genotype profile in JSON format")
    cli.add_argument("profile2", help="simulated or inferred genotype profile in JSON format")


def main(args):
    contained, total = mhapi.contain(
        Profile(fromfile=args.profile1), Profile(fromfile=args.profile2)
    )
    data = {
        "containment": round(contained / total, 4),
        "contained_alleles": contained,
        "total_alleles": total,
    }
    with mhopen(args.out, "w") as fh:
        json.dump(data, fh, indent=4)
