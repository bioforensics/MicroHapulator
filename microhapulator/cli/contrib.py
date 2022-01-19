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
from microhapulator.parsers import open as mhopen
import microhapulator.api as mhapi
from microhapulator.profile import Profile


def subparser(subparsers):
    desc = "Estimate the minimum number of DNA contributors to a suspected mixture"
    cli = subparsers.add_parser("contrib", description=desc)
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument("result", help="MicroHapulator typing result in JSON format")


def main(args):
    ncontrib, nloci, ploci = mhapi.contrib(Profile(fromfile=args.result))
    data = {
        "min_num_contrib": ncontrib,
        "num_loci_max_alleles": nloci,
        "perc_loci_max_alleles": ploci,
    }
    with mhopen(args.out, "w") as fh:
        json.dump(data, fh, indent=4)
