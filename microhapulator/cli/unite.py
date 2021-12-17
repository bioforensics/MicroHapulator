#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser("unite")
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "-s",
        "--seed",
        type=int,
        default=None,
        metavar="INT",
        help="seed for " "random number generator",
    )
    cli.add_argument("mom", help="simulated or inferred genotype in JSON format")
    cli.add_argument("dad", help="simulated or inferred genotype in JSON format")
