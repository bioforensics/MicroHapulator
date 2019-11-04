#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    desc = (
        'Use precomputed population allele frequencies to simulate a diploid '
        'DNA profile for the specified panel of microhaplotypes.'
    )
    cli = subparsers.add_parser('sim', description=desc)

    cli.add_argument(
        '-o', '--out', metavar='FILE', help='write simulated profile data in '
        'JSON format to FILE'
    )
    cli.add_argument(
        '--haplo-seq', metavar='FILE', help='write simulated haplotype '
        'sequences in FASTA format to FILE'
    )
    cli.add_argument(
        '-s', '--seed', type=int, default=None, metavar='INT', help='seed for '
        'random number generator'
    )
    cli.add_argument(
        '-r', '--relaxed', action='store_true', help='if a marker in the panel '
        'has no frequency data for a requested population, randomly draw an '
        'allele (from a uniform distribution) from all possible alleles; by '
        'default, these markers are exluded from simulation'
    )
    cli.add_argument(
        'popid', nargs=2, help='population identifiers for the two parental '
        'haplotypes'
    )
    cli.add_argument(
        'panel', nargs='+', help='panel from which to simulate microhap '
        'profiles; can be a list of MicroHapDB marker identifiers or a '
        'filename containing marker identifiers (one per line)'
    )
