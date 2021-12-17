# -----------------------------------------------------------------------------
# Copyright (c) 2021, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import sys


def subparser(subparsers):
    cli = subparsers.add_parser("balance")
    cli.add_argument("-c", "--csv", metavar="FILE", help="write read counts to FILE in CSV format")
    cli.add_argument(
        "-D",
        "--no-discarded",
        dest="discarded",
        action="store_false",
        help="do not included mapping but discarded reads in read counts; by default, reads that "
        "are mapped to the marker but discarded because they do not span all variants at the "
        "marker are included",
    )
    cli.add_argument("input", help="typing result in JSON format")
