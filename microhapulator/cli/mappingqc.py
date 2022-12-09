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
        "marker", help="path of csv file containing number of reads mapped to marker sequences"
    )
    cli.add_argument(
        "refr", help="path of csv file containing number of reads mapped to full reference genome"
    )
    cli.add_argument(
        "rep", help="path of csv file containing number of repetitive reads per marker"
    )
    cli.add_argument("csv", help="write read counts to FILE in CSV format")
    cli.add_argument(
        "figure",
        help="create donut plot to FILE showing porportions of on target, off target, repetitive, and contaminant reads ",
    )


def main(args):
    data = mhapi.read_mapping_qc(args.marker, args.refr, args.rep, args.figure)
    data.to_csv(args.csv, index=False)
