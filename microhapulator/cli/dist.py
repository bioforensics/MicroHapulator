#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser('dist')
    cli.add_argument(
        '-o', '--out', metavar='FILE', help='write output to "FILE"; by '
        'default, output is written to the terminal (standard output)'
    )
    cli.add_argument(
        '-p', '--ploidy', type=int, metavar='PLD', default=2,
        help='for simulated genotypes, the number of distinct haplotypes '
        'present in the sample; default is 2 (for a single diploid genotype)'
    )
    cli.add_argument(
        '--sim', nargs='*', default=[], help='simulated genotype data in BED format'
    )
    cli.add_argument(
        '--obs', nargs='*', default=[], help='observed genotype data in JSON format'
    )
