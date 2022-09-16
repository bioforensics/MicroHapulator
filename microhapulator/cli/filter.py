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
from microhapulator.parsers import load_marker_thresholds
from microhapulator.profile import TypingResult
import pandas as pd
import sys


def subparser(subparsers):
    desc = "Apply static and/or dynamic thresholds to distinguish true and false haplotypes. Thresholds are applied to the haplotype read counts of a raw typing result. Static integer thresholds are commonly used as detection thresholds, below which any haplotype count is considered noise. Dynamic thresholds are commonly used as analytical thresholds and represent a percentage of the total read count at the marker, after any haplotypes failing a static threshold are discarded."
    cli = subparsers.add_parser("filter", description=desc)
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        default=sys.stdout,
        help="write output to FILE; by default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "-s",
        "--static",
        metavar="ST",
        type=int,
        default=5,
        help="global fixed read count threshold",
    )
    cli.add_argument(
        "-d",
        "--dynamic",
        metavar="DT",
        type=float,
        default=0.02,
        help="global percentage of total read count threshold; e.g. use --dynamic=0.02 to apply a 2%% analytical threshold",
    )
    cli.add_argument(
        "-c",
        "--config",
        metavar="CSV",
        default=None,
        help="CSV file specifying marker-specific thresholds to override global thresholds; three required columns: 'Marker' for the marker name; 'Static' and 'Dynamic' for marker-specific thresholds",
    )
    cli.add_argument("result", help="MicroHapulator typing result in JSON format")


def main(args):
    result = TypingResult(fromfile=args.result)
    thresholds = load_marker_thresholds(result.markers(), args.config, args.static, args.dynamic)
    result.filter(thresholds)
    result.dump(args.out)
