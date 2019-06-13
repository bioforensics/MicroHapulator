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
        'genotype for the specified panel of microhaplotypes.'
    )
    cli = subparsers.add_parser('sim', description=desc)

    cli.add_argument(
        '-o', '--out', metavar='FILE', help='write simulated genotype data in '
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
        '-r', '--relaxed', action='store_true', help='if a locus in the panel '
        'has no frequency data for a requested population, randomly draw an '
        'allele (from a uniform distribution) from all possible alleles; by '
        'default, these loci are exluded from simulation'
    )
    cli.add_argument(
        'popid', nargs=2, help='population identifiers for two haplotypes'
    )
    cli.add_argument(
        'panel', nargs='+', help='panel from which to simulate microhap '
        'genotypes; can be a list of MicroHapDB locus identifiers or the '
        'label for one of MicroHapulator\'s preset panels (such as alpha, '
        'beta, or usa)'
    )
