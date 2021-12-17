#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser("contrib")
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "-j",
        "--json",
        metavar="FILE",
        help="precomputed genotype profile in "
        "JSON format; if not specified, must supply arguments for "
        "`--refr-fasta` and `--bam` flags",
    )
    cli.add_argument(
        "-r", "--refr", metavar="FILE", help="microhap marker sequences in " "Fasta format"
    )
    cli.add_argument(
        "-b", "--bam", metavar="FILE", help="aligned and sorted reads in BAM " "format"
    )
    cli.add_argument(
        "-s",
        "--static",
        metavar="ST",
        type=int,
        default=None,
        help="apply a static threshold for calling genotypes; see `mhpl8r type --help`",
    )
    cli.add_argument(
        "-d",
        "--dynamic",
        metavar="DT",
        type=float,
        default=None,
        help="apply a dynamic threshold for calling genotypes; see `mhpl8r type --help`",
    )
