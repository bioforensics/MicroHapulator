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
    desc = "Calculate number of on target, off target, repetitive, and contaminant reads and create a donut plot"
    cli = subparsers.add_parser("mappingqc", description=desc)
    cli.add_argument(
        "--marker",
        help="path of csv file containing number of reads mapped to marker sequences",
        required=True,
    )
    cli.add_argument(
        "--refr",
        help="path of csv file containing number of reads mapped to full reference genome",
        required=True,
    )
    cli.add_argument(
        "--rep",
        help="path of csv file containing number of repetitive reads per marker",
        required=True,
    )
    cli.add_argument("--csv", help="write read counts to FILE in CSV format", required=True)
    cli.add_argument(
        "--figure",
        help="create donut plot to FILE showing porportions of on target, off target, repetitive, and contaminant reads",
        required=True,
    )
    cli.add_argument(
        "--title",
        default=None,
        help="add a title (such as a sample name) to the histogram plot",
    )


def main(args):
    data = mhapi.read_mapping_qc(args.marker, args.refr, args.rep, args.figure, args.title)
    data.to_csv(args.csv, index=False)
