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
import sys


def subparser(subparsers):
    desc = "Calculate off target read mapping"
    cli = subparsers.add_parser("offtarget", description=desc)
    cli.add_argument("markerbam", help="alignment file of reads aligned to marker sequences")
    cli.add_argument(
        "refbam", help="alignment file in BAM format of reads aligned to complete reference genome"
    )
    cli.add_argument(
        "tsv",
        help="marker definitions tsv including chromosome and full reference genome offset columns",
    )
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        default=sys.stdout,
        help="write output to FILE; by default, output is written to the terminal (standard output)",
    )


def main(args):
    data = mhapi.off_target_mapping(args.markerbam, args.refbam, args.tsv)
    data.to_csv(args.out, index=False)
