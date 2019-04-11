#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser('refr')
    cli.add_argument(
        '-o', '--out', metavar='FILE', default='-',
        help='write output to "FILE"; by default, output is written to the '
        'terminal (standard output)'
    )
    cli.add_argument(
        '-d', '--delta', type=int, default=30, metavar='Δ',
        help='extend each microhap locus by Δ nucleotides'
    )
    cli.add_argument(
        '-m', '--min-length', type=int, default=350, metavar='M',
        help='after applying deltas, if a microhap locus is shorter than M '
        'nucleotides, extend both sides equally so that it is M nucleotides '
        'in length'
    )
    cli.add_argument(
        'refrfasta', help='reference genome file'
    )
    cli.add_argument(
        'panel', nargs='*', help='list of MicroHapDB locus IDs; by default, a '
        'panel of 22 ALFRED microhaplotype loci is used'
    )
