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
import microhapulator.api as mhapi
from microhapulator.profile import Profile
from pkg_resources import resource_filename
from snakemake import snakemake


def subparser(subparsers):
    desc = "Perform a complete end-to-end microhap analysis pipeline"
    cli = subparsers.add_parser("pipe", description=desc)
    cli.add_argument(
        "-w",
        "--workdir",
        metavar="D",
        help="pipeline working directory; default is current directory",
    )
    cli.add_argument("seqpath", help="path to a directory containing FASTQ files")
    cli.add_argument("samples", help="sample names")


def main(args):
    config = dict(
        samples=sorted(args.samples),
        readfiles=["README.md"],  # FIXME,
        mhrefr="README.md",
        mhdefn="README.md",
    )
    snakefile = resource_filename("microhapulator", "Snakefile")
    snakemake(snakefile, dryrun=True, config=config)
