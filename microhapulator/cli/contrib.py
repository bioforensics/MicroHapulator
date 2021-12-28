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
from microhapulator import open as mhopen
from microhapulator.op.contrib import contrib
from microhapulator.profile import Profile
from microhapulator.op.type import type as mhtype


def subparser(subparsers):
    cli = subparsers.add_parser("contrib")
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument(
        "-j",
        "--json",
        metavar="FILE",
        help="precomputed genotype profile in "
        "JSON format; if not specified, must supply arguments for "
        "`--refr-fasta` and `--bam` flags",
    )
    cli.add_argument(
        "-t", "--tsv", metavar="FILE", help="microhap marker definitions in tabular (TSV) format"
    )
    cli.add_argument("-b", "--bam", metavar="FILE", help="aligned and sorted reads in BAM format")
    cli.add_argument(
        "-s",
        "--static",
        metavar="ST",
        type=int,
        default=None,
        help="apply a static threshold for calling genotypes; see `mhpl8r type --help`",
    )
    cli.add_argument(
        "-d",
        "--dynamic",
        metavar="DT",
        type=float,
        default=None,
        help="apply a dynamic threshold for calling genotypes; see `mhpl8r type --help`",
    )


def load_profile(bamfile=None, markertsv=None, json=None, **kwargs):
    if not json and (not bamfile or not markertsv):
        message = "must provide either JSON profile or BAM and refr FASTA"
        raise ValueError(message)
    if json:
        profile = Profile(fromfile=json)
    else:
        profile = mhtype(bamfile, markertsv, **kwargs)
    return profile


def main(args):
    profile = load_profile(
        bamfile=args.bam,
        markertsv=args.tsv,
        json=args.json,
        static=args.static,
        dynamic=args.dynamic,
    )
    ncontrib, nloci, ploci = contrib(profile)
    data = {
        "min_num_contrib": ncontrib,
        "num_loci_max_alleles": nloci,
        "perc_loci_max_alleles": ploci,
    }
    with mhopen(args.out, "w") as fh:
        json.dump(data, fh, indent=4)
