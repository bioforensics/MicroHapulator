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
        'If a single profile is provided, the random match probability of '
        'the profile is computed. If a pair of profiles is provided, a '
        'likelihood ratio test is performed comparing the likelihood that the '
        'two profiles are from the same individual versus the likelihood '
        'that the two profiles are from random unrelated individuals. The '
        'genotype profiles are assumed to be identical, and differences '
        'between the two profiles are assumed to be the result of genotyping '
        'error. The test does not make sense for profiles with many allele '
        'differences.'
    )
    cli = subparsers.add_parser('prob', description=desc)
    cli.add_argument(
        '-e', '--erate', type=float, metavar='ε', default=0.001, help='rate '
        'at which errors in genotyping are expected; default is 0.001'
    )
    cli.add_argument(
        '-o', '--out', metavar='FILE', help='write output to "FILE"; by '
        'default, output is written to the terminal (standard output)'
    )
    cli.add_argument(
        'population', help='indicate which allele frequencies to use'
    )
    cli.add_argument(
        'profile1', help='profile in JSON format'
    )
    cli.add_argument(
        'profile2', nargs='?', default=None, help='profile in JSON format'
    )
