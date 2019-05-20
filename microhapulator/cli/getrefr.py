#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser('getrefr')
    cli.add_argument('--debug', action='store_true', help='debugging output')
    cli.add_argument('filepath', nargs='?', help='path to a local copy of the '
                     'GRCh38 reference genome; if not provided, it is '
                     'downloaded from UCSC (requires an internet connection)')
