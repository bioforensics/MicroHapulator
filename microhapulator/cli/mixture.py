#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser('mixture')
    cli.add_argument(
        '--indiv', nargs='+', metavar='POP', action='append',
        help='population ID (or pair of IDs) for each individual; invoke the '
        '`--indiv` flag once for each individual in the mixture'
    )
    cli.add_argument(
        '--panel', nargs='+', metavar='ID', help='list of MicroHapDB locus '
        'IDs for which to simulate data; by default, a panel of 22 ALFRED '
        'microhaplotype loci is used'
    )
    cli.add_argument(
        '-r', '--relaxed', action='store_true', help='if a locus in the panel '
        'has no frequency data for a requested population(s), randomly draw '
        'an allele (from a uniform distribution) from all possible alleles; '
        'by default, these loci are exluded from simulation'
    )
    cli.add_argument(
        '-p', '--proportions', type=float, nargs='+', metavar='P',
        help='simulate mixture contributors at the specified proportions; by '
        'default even proportions are used'
    )
    cli.add_argument(
        '--hap-seeds', type=int, default=None, metavar='INT', nargs='+',
        help='random seeds for simulating haplotypes'
    )
    cli.add_argument(
        '--genotype', metavar='FILE', help='write simulated genotype data in '
        'JSON format to FILE'
    )
    cli.add_argument(
        '-n', '--num-reads', type=int, default=500000, metavar='N',
        help='number of reads to simulate; default is 500000'
    )
    cli.add_argument(
        '--seq-seeds', type=int, default=None, metavar='INT', nargs='+',
        help='random seeds for simulated sequencing'
    )
    cli.add_argument(
        '--seq-threads', type=int, default=None, metavar='INT', help='number '
        'of threads to use when simulating targeted amplicon sequencing'
    )
    cli.add_argument(
        '-y', '--dry-run', action='store_true', help='simulate genotype and '
        'produce any requested outputs, but terminate prior to simulating '
        'reads'
    )
    cli.add_argument(
        '-o', '--out', metavar='FILE', help='write simulated MiSeq reads in '
        'FASTQ format to FILE; by default, reads are written to terminal '
        '(standard output)'
    )
