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


import microhapulator.api as mhapi
from microhapulator.profile import TypingResult


def subparser(subparsers):
    desc = (
        "Plot interlocus balance in the terminal and/or a high-resolution graphic. Also normalize "
        "read counts and perform a chi-square goodness-of-fit test assuming uniform read coverage "
        "across markers. The reported chi-square statistic measures the extent of imbalance, and "
        "can be compared among samples sequenced using the same panel: the minimum value of 0 "
        "represents perfectly uniform coverage, while the maximum value of D occurs when all "
        "reads map to a single marker (D represents the degrees of freedom, or the number of "
        "markers minus 1)."
    )
    cli = subparsers.add_parser("balance", description=desc)
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
    cli.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="do not print interlocus balance histogram to standard output in ASCII",
    )
    cli.add_argument(
        "--figure",
        metavar="FILE",
        default=None,
        help="plot interlocus balance histogram to FILE using Matplotlib; image format is inferred from extension of provided file name",
    )
    cli.add_argument(
        "--figsize",
        metavar=("W", "H"),
        nargs=2,
        type=float,
        default=(6, 4),
        help="dimensions (width × height in inches) of the image file to be generated; 6 4 by default",
    )
    cli.add_argument(
        "--dpi",
        metavar="DPI",
        type=int,
        default=200,
        help="resolution (in dots per inch) of the image file to be generated; DPI=200 by default",
    )
    cli.add_argument(
        "--color",
        metavar="COL",
        default="#1f77b4",
        help="color of the histogram to be generated in the image file; COL='#1f77b4' by default",
    )
    cli.add_argument("input", help="a typing result including haplotype counts in JSON format")


def main(args):
    result = TypingResult(fromfile=args.input)
    chisq, data = mhapi.balance(
        result,
        include_discarded=args.discarded,
        terminal=not args.quiet,
        tofile=args.figure,
        figsize=args.figsize,
        dpi=args.dpi,
        color=args.color,
    )
    print(f"Extent of imbalance (chi-square statistic): {chisq:.4f}")
    if args.csv:
        data.to_csv(args.csv, index=False)
