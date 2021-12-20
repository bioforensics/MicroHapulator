#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    desc = "Simulate a diploid genotype using the given population microhaplotype frequencies"
    cli = subparsers.add_parser("sim", description=desc)
    cli.add_argument(
        "freq", help="population microhaplotype frequencies in tabular (tab separated) format"
    )
    cli.add_argument(
        "-s",
        "--seed",
        type=int,
        default=None,
        metavar="INT",
        help="seed for random number generator",
    )
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help="write simulated profile data in JSON format to FILE",
    )
    cli.add_argument(
        "--haplo-seq",
        metavar="FILE",
        help="write simulated haplotype sequences in FASTA format to FILE",
    )
    cli.add_argument(
        "--sequences",
        metavar="FILE",
        help="microhaplotype sequences in FASTA format; required if `--haplo-seq` enabled, ignored if not",
    )
    cli.add_argument(
        "--markers",
        metavar="FILE",
        help="microhaplotype marker definitions in tabular (tab separated) format; required if `--haplo-seq` enabled, ignored if not",
    )
