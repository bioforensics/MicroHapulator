# -------------------------------------------------------------------------------------------------
# Copyright (c) 2021, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from microhapulator.op import balance
from microhapulator.profile import ObservedProfile


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


def main(args):
    profile = ObservedProfile(fromfile=args.input)
    data = balance(profile, include_discarded=args.discarded)
    if args.csv:
        data.to_csv(args.csv, index=False)
