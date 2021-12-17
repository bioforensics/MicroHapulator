#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import SUPPRESS
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
        "-e",
        "--effcov",
        metavar="EC",
        type=float,
        default=0.25,
        help="only reads that span all SNPs in a microhaplotype are retained, all others are "
        "discarded; if most of the reads related to a marker are discarded, it has low *effective "
        "coverage*; this parameter sets that threshold, i.e., if the fraction of retained reads "
        "is < EC it is considered low effective coverage; by default EC=0.25 (or 25%%)",
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
        "defined, --static is only applied to markers with low effective coverage",
    )
    cli.add_argument(
        "-d",
        "--dynamic",
        metavar="DT",
        type=float,
        default=None,
        help="apply a dynamic threshold for calling genotypes, i.e., if AC is the average count "
        "of all haplotypes observed at the marker, discard any haplotypes whose count is less "
        "than DT * AC; by default, DT is undefined, the dynamic filter is not applied, and only "
        "raw haplotype counts are reported, not genotype calls; if --static is also defined, "
        "--dynamic is only applied to markers with high effective coverage",
    )
    cli.add_argument("-m", "--max-depth", metavar="M", type=float, default=1e6, help=SUPPRESS)
    cli.add_argument("tsv", help="microhap marker definitions in TSV format")
    cli.add_argument("bam", help="aligned and sorted reads in BAM format")
