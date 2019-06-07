#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


def subparser(subparsers):
    cli = subparsers.add_parser('sim')
    cli._positionals.title = 'Input configuration'
    cli._optionals.title = 'Miscellaneous'

    hapargs = cli.add_argument_group('Haplotype simulation')
    hapargs.add_argument(
        '--panel', nargs='+', metavar='ID', help='list of MicroHapDB locus '
        'IDs for which to simulate data; by default, a panel of 22 ALFRED '
        'microhaplotype loci is used'
    )
    hapargs.add_argument(
        '-r', '--relaxed', action='store_true', help='if a locus in the panel '
        'has no frequency data for a requested population, randomly draw an '
        'allele (from a uniform distribution) from all possible alleles; by '
        'default, these loci are exluded from simulation'
    )
    hapargs.add_argument(
        '--hap-seed', type=int, default=None, metavar='INT', help='random '
        'seed for simulating haplotypes'
    )

    seqargs = cli.add_argument_group('Targeted sequencing')
    seqargs.add_argument(
        '-n', '--num-reads', type=int, default=500000, metavar='N',
        help='number of reads to simulate; default is 500000'
    )
    seqargs.add_argument(
        '--seq-seed', type=int, default=None, metavar='INT', help='random '
        'seed for simulated sequencing'
    )
    seqargs.add_argument(
        '--seq-threads', type=int, default=None, metavar='INT', help='number '
        'of threads to use when simulating targeted amplicon sequencing'
    )
    outargs = cli.add_argument_group('Output configuration')
    outargs.add_argument(
        '-y', '--dry-run', action='store_true', help='simulate genotype and '
        'produce any requested outputs, but terminate prior to simulating '
        'reads'
    )
    outargs.add_argument(
        '-o', '--out', metavar='FILE', help='write simulated MiSeq reads in '
        'FASTQ format to FILE; by default, reads are written to terminal '
        '(standard output)'
    )
    outargs.add_argument(
        '--genotype', metavar='FILE', help='write simulated genotype data in '
        'JSON format to FILE'
    )
    outargs.add_argument(
        '--haploseq', metavar='FILE', help='write simulated haplotype '
        'sequences in FASTA format to FILE'
    )

    cli.add_argument('popid', nargs='+', help='population ID(s)')
    cli._action_groups[1], cli._action_groups[-1] = cli._action_groups[-1], cli._action_groups[1]
