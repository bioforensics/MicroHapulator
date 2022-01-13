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


import microhapulator.api as mhapi
import sys


def subparser(subparsers):
    cli = subparsers.add_parser("type")
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        default=sys.stdout,
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "-b",
        "--base-qual",
        metavar="B",
        type=int,
        default=10,
        help="minimum base quality required for haplotype calling; by default B=10, "
        "corresponding to Q10, i.e., 90%% probability that base call is correct",
    )
    cli.add_argument(
        "-m",
        "--max-depth",
        metavar="M",
        type=float,
        default=1e6,
        help="maximum permitted read depth; by default M=1000000",
    )
    cli.add_argument("tsv", help="microhap marker definitions in TSV format")
    cli.add_argument("bam", help="aligned and sorted reads in BAM format")


def main(args):
    profile = mhapi.type(args.bam, args.tsv, minbasequal=args.base_qual, max_depth=args.max_depth)
    profile.dump(args.out)
