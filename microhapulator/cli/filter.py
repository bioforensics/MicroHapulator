# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


import microhapulator.api as mhapi
from microhapulator.profile import TypingResult
import pandas as pd
import sys


def subparser(subparsers):
    cli = subparsers.add_parser("filter")
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        default=sys.stdout,
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "-s",
        "--static",
        metavar="ST",
        type=int,
        default=None,
        help="apply a static threshold for calling genotypes, i.e., discard any haplotype whose "
        "count is less than ST; by default, ST is undefined, the static filter is not applied, "
        "and only raw haplotype counts are reported, not genotype calls; if --dynamic is also "
        "defined, --static is applied first",
    )
    cli.add_argument(
        "-d",
        "--dynamic",
        metavar="DT",
        type=float,
        default=None,
        help="apply a dynamic threshold for calling genotypes, i.e., if C is the total count "
        "of all haplotypes observed at the marker, discard any haplotypes whose count is less "
        "than C * DT; by default, DT is undefined, the dynamic filter is not applied, and only "
        "raw haplotype counts are reported, not genotype calls; if --static is also defined, "
        "it is applied first, and the counts of any haplotypes that do not pass the --static "
        "filter are deducted from the total count C",
    )
    cli.add_argument(
        "-c",
        "--config",
        metavar="FILE",
        default=None,
        help="CSV file with marker-specific static and dynamic thresholds; the file should "
        "contain 3 columns named Marker, Static, and Dynamic; there should be at most one row per "
        "marker; if --static and/or --dynamic are specified, these serve as default values, while "
        "any thresholds defined in the config file will override these values for the specified "
        "markers",
    )
    cli.add_argument("result", help="MicroHapulator typing result in JSON format")


def main(args):
    result = TypingResult(fromfile=args.result)
    config = None
    if args.config:
        config = pd.read_csv(args.config, sep=None, engine="python")
    result.filter(static=args.static, dynamic=args.dynamic, config=config)
    result.dump(args.out)
