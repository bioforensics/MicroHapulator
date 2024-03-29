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


from microhapulator import open as mhopen
from microhapulator.profile import SimulatedProfile


def subparser(subparsers):
    desc = "Combine simulated profiles into a mock DNA mixture"
    cli = subparsers.add_parser("mix", description=desc)
    cli.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help='write output to "FILE"; by '
        "default, output is written to the terminal (standard output)",
    )
    cli.add_argument("profiles", nargs="+", help="simulated genotype profiles in JSON format")


def main(args):
    profiles = [SimulatedProfile(pfile) for pfile in args.profiles]
    combined = SimulatedProfile.merge(profiles)
    with mhopen(args.out, "w") as fh:
        combined.dump(fh)
