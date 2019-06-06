#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser('type')
    cli.add_argument(
        '-o', '--out', metavar='FILE', help='write output to "FILE"; by '
        'default, output is written to the terminal (standard output)'
    )
    cli.add_argument(
        '-t', '--threshold', type=int, default=None, help='coverage threshold '
        'below which alleles will be ignored; by default, threshold is '
        'computed automatically'
    )
    cli.add_argument(
        'refr', help='microhap locus sequences in Fasta format'
    )
    cli.add_argument(
        'bam', help='aligned and sorted reads in BAM format'
    )
