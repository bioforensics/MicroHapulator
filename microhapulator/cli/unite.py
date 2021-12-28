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

import numpy.random
from microhapulator import open as mhopen
from microhapulator.profile import Profile


def subparser(subparsers):
    cli = subparsers.add_parser("unite")
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "-s",
        "--seed",
        type=int,
        default=None,
        metavar="INT",
        help="seed for " "random number generator",
    )
    cli.add_argument("mom", help="simulated or inferred genotype in JSON format")
    cli.add_argument("dad", help="simulated or inferred genotype in JSON format")


def main(args):
    if args.seed:
        numpy.random.seed(args.seed)
    profile = Profile.unite(
        Profile(fromfile=args.mom),
        Profile(fromfile=args.dad),
    )
    with mhopen(args.out, "w") as fh:
        profile.dump(fh)
