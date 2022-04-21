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


def subparser(subparsers):
    desc = "Compute and plot heterozygote balance"
    cli = subparsers.add_parser("hetbalance", description=desc)
    cli.add_argument("-c", "--csv", metavar="FILE", help="write read counts to FILE in CSV format")
    cli.add_argument(
        "--figure",
        metavar="FILE",
        default=None,
        help="plot heterzygote balance bar graph to FILE using Matplotlib; image format is inferred from extension of provided file name",
    )
    cli.add_argument(
        "--figsize",
        metavar=("W", "H"),
        nargs=2,
        type=float,
        default=None,
        help="dimensions (width × height in inches) of the image file to be generated; figure dimensions determined automatically by default",
    )
    cli.add_argument(
        "--dpi",
        metavar="DPI",
        type=int,
        default=200,
        help="resolution (in dots per inch) of the image file to be generated; DPI=200 by default",
    )
    cli.add_argument(
        "-t",
        "--title",
        metavar="T",
        default=None,
        help="add a title (such as a sample name) to the histogram plot",
    )
    cli.add_argument(
        "--labels", action="store_true", help="include labels showing marker names and read counts"
    )
    cli.add_argument(
        "--absolute", action="store_true", help="plot absolute rather than relative read counts"
    )
    cli.add_argument("input", help="a typing result including haplotype counts in JSON format")


def main(args):
    result = TypingResult(fromfile=args.input)
    tstat, data = mhapi.heterozygote_balance(
        result,
        tofile=args.figure,
        title=args.title,
        figsize=args.figsize,
        dpi=args.dpi,
        dolabels=args.labels,
        absolute=args.absolute,
    )
    print(f"Extent of imbalance (t-statistic): {tstat:.4f}")
    if args.csv:
        data.to_csv(args.csv, index=False)
