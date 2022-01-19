# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


import microhapulator.api as mhapi
from microhapulator.profile import TypingResult
import pandas as pd
import sys


def subparser(subparsers):
    desc = "Convert a typing result to a format compatible with probabilistic genotyping software applications"
    cli = subparsers.add_parser("convert", description=desc)
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        default=sys.stdout,
        help="write output to 'FILE'; by default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "--no-counts",
        dest="counts",
        action="store_false",
        help="do not include haplotype counts if you are interpreting your data with a semi-"
        "continuous probgen model such as LRMix Studio; by default, haplotype counts are included "
        "for interpretation with fully continuous probgen model such as EuroForMix",
    )
    cli.add_argument(
        "-f",
        "--fix-homo",
        action="store_true",
        help="duplicate a homozygous haplotype so that it is reported twice",
    )
    cli.add_argument("result", help="filtered MicroHapulator typing result in JSON format")
    cli.add_argument("sample", help="sample name")


def main(args):
    result = TypingResult(fromfile=args.result)
    result.dump_csv(args.out, args.sample, counts=args.counts, fix_homo=args.fix_homo)
