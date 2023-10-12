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
    cli = subparsers.add_parser("type", description="Perform haplotype calling")
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        default=sys.stdout,
        help="write output to FILE; by default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "-b",
        "--base-qual",
        metavar="B",
        type=int,
        default=10,
        help="minimum base quality (PHRED score) to be considered reliable for haplotype calling; by default B=10, corresponding to Q10, i.e., 90%% probability that the base call is correct",
    )
    cli.add_argument(
        "-m",
        "--max-depth",
        metavar="M",
        type=float,
        default=1e6,
        help="maximum permitted read depth; by default M=1000000",
    )
    cli.add_argument(
        "-r",
        "--relaxed",
        dest="strict",
        action="store_false",
        help="disable strict validation of marker definitions",
    )
    cli.add_argument(
        "tsv",
        help="path of a TSV file containing marker metadata, specifically the offset of each SNP for every marker in the panel",
    )
    cli.add_argument(
        "bam",
        help="path of a BAM file containing NGS reads aligned to marker reference sequences and sorted",
    )


def main(args):
    profile = mhapi.type(
        args.bam,
        args.tsv,
        minbasequal=args.base_qual,
        max_depth=args.max_depth,
        strict=args.strict,
    )
    profile.dump(args.out)
