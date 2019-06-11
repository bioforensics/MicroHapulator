#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse


def subparser(subparsers):
    desc = (
        'Given one or more simulated diploid genotypes, simulate Illumina '
        'MiSeq sequencing of the given sample.'
    )
    cli = subparsers.add_parser('seq', description=desc)

    cli.add_argument(
        '-o', '--out', metavar='FILE', help='write simulated MiSeq reads in '
        'FASTQ format to FILE; by default, reads are written to terminal '
        '(standard output)'
    )
    cli.add_argument(
        '-n', '--num-reads', type=int, default=500000, metavar='N',
        help='number of reads to simulate; default is 500000'
    )
    cli.add_argument(
        '--threads', type=int, default=None, metavar='INT', help='number of '
        'threads to use when simulating targeted amplicon sequencing'
    )
    cli.add_argument(
        '-p', '--proportions', type=float, nargs='+', metavar='P',
        help='simulated mixture samples with multiple contributors at the '
        'specified proportions; by default even proportions are used'
    )
    cli.add_argument(
        '-s', '--seeds', nargs='+', type=int, default=None, metavar='INT',
        help='seeds for random number generator, 1 per genotype'
    )
    cli.add_argument('--signature', default=None, help=argparse.SUPPRESS)
    cli.add_argument(
        '-d', '--debug', action='store_true', help='print debugging info'
    )
    cli.add_argument(
        'genotypes', nargs='+', help='one or more simple or complex genotype '
        'JSON files'
    )
