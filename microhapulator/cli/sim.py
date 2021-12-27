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


def subparser(subparsers):
    desc = "Simulate a diploid genotype using the given population microhaplotype frequencies"
    cli = subparsers.add_parser("sim", description=desc)
    cli.add_argument(
        "freq", help="population microhaplotype frequencies in tabular (tab separated) format"
    )
    cli.add_argument(
        "-s",
        "--seed",
        type=int,
        default=None,
        metavar="INT",
        help="seed for random number generator",
    )
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help="write simulated profile data in JSON format to FILE",
    )
    cli.add_argument(
        "--haplo-seq",
        metavar="FILE",
        help="write simulated haplotype sequences in FASTA format to FILE",
    )
    cli.add_argument(
        "--sequences",
        metavar="FILE",
        help="microhaplotype sequences in FASTA format; required if `--haplo-seq` enabled, ignored if not",
    )
    cli.add_argument(
        "--markers",
        metavar="FILE",
        help="microhaplotype marker definitions in tabular (tab separated) format; required if `--haplo-seq` enabled, ignored if not",
    )
