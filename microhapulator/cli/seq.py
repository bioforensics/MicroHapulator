#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import Action, SUPPRESS
import microhapulator
import sys


class OutfilesAction(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        numargs = len(values)
        if numargs < 1 or numargs > 2:
            message = f'expected 1 or 2 output filenames, got {numargs}: {",".join(values)}'
            parser.error(message)
        filehandles = [open(filename, "w") for filename in values]
        setattr(namespace, self.dest, filehandles)


def subparser(subparsers):
    desc = (
        "Given one or more diploid genotype profiles, simulate Illumina MiSeq "
        'sequencing of the given "sample."'
    )
    cli = subparsers.add_parser("seq", description=desc)
    cli.add_argument(
        "-o",
        "--out",
        nargs="+",
        default=[sys.stdout],
        action=OutfilesAction,
        help="write simulated paired-end MiSeq reads in FASTQ format to the specified file(s); if "
        "one filename is provided, reads are interleaved and written to the file; if two "
        "filenames are provided, reads are written to paired files; by default, reads are "
        "interleaved and written to the terminal (standard output)",
    )
    cli.add_argument(
        "-n",
        "--num-reads",
        type=int,
        default=500000,
        metavar="N",
        help="number of reads to simulate; default is 500000",
    )
    cli.add_argument(
        "-p",
        "--proportions",
        type=float,
        nargs="+",
        metavar="P",
        help="simulated mixture samples with multiple contributors at the "
        "specified proportions; by default even proportions are used",
    )
    cli.add_argument(
        "-s",
        "--seeds",
        nargs="+",
        type=int,
        default=None,
        metavar="INT",
        help="seeds for random number generator, 1 per profile",
    )
    cli.add_argument("--signature", default=None, help=SUPPRESS)
    cli.add_argument("--threads", type=int, default=1, metavar="INT", help=SUPPRESS)
    cli.add_argument(
        "profiles", nargs="+", help="one or more simple or complex profiles " "(JSON files)"
    )
