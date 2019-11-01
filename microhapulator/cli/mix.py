#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser('mix')
    cli.add_argument(
        '-o', '--out', metavar='FILE', help='write output to "FILE"; by '
        'default, output is written to the terminal (standard output)'
    )
    cli.add_argument(
        'profiles', nargs='+', help='simulated genotype profiles in JSON format'
    )
