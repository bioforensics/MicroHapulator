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
    desc = "Calculate number of reads that map to a marker sequence but map preferentially to another locus when aligned to the whole genome"
    cli = subparsers.add_parser("repetitive", description=desc)
    cli.add_argument("markerbam", help="alignment file of reads aligned to marker sequences")
    cli.add_argument("refbam", help="alignment file in BAM format of reads aligned to hg38")
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
    cli.add_argument(
        "-b",
        "--base-qual",
        metavar="B",
        type=int,
        default=10,
        help="minimum base quality (PHRED score) to be considered reliable for haplotype calling; by default B=10, corresponding to Q10, i.e., 90%% probability that the base call is correct",
    )


def main(args):
    data = mhapi.repetitive_mapping(
        args.markerbam, args.refbam, args.tsv, minbasequal=args.base_qual
    )
    data.to_csv(args.out, index=False)
