# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


import json
from microhapulator.parsers import load_marker_frequencies
from microhapulator.parsers import open as mhopen
import microhapulator.api as mhapi
from microhapulator.profile import Profile


def subparser(subparsers):
    desc = "Compute a profile random match probability (RMP) or an RMP-based likelihood ratio (LR) test"
    cli = subparsers.add_parser("prob", description=desc)
    cli.add_argument(
        "-e",
        "--erate",
        type=float,
        metavar="ε",
        default=0.001,
        help="rate of genotyping error; by default ε=0.001",
    )
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help="write output to FILE; by default, output is written to the terminal (standard output)",
    )
    cli.add_argument("freq", help="population haplotype frequencies in tabular (TSV) format")
    cli.add_argument("profile1", help="typing result or simulated genotype in JSON format")
    cli.add_argument(
        "profile2",
        nargs="?",
        default=None,
        help="typing result or simulated genotype in JSON format; optional",
    )


def main(args):
    prof1 = Profile(fromfile=args.profile1)
    prof2 = Profile(fromfile=args.profile2) if args.profile2 else None
    frequencies = load_marker_frequencies(args.freq)
    result = mhapi.prob(frequencies, prof1, prof2=prof2, erate=args.erate)
    key = "random_match_probability" if prof2 is None else "likelihood_ratio"
    data = {
        key: "{:.3E}".format(result),
    }
    with mhopen(args.out, "w") as fh:
        json.dump(data, fh, indent=4)
