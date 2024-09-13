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


from microhapulator import load_marker_frequencies
import microhapulator.api as mhapi
from microhapulator.marker import MicrohapIndex
import sys


def subparser(subparsers):
    desc = "Simulate a diploid genotype from the specified microhaplotype frequencies"
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


def main(args):
    frequencies = load_marker_frequencies(args.freq)
    profile = mhapi.sim(frequencies, seed=args.seed)
    with open(args.out, "w") as fh:
        profile.dump(fh)
        message = "profile JSON written to {:s}".format(fh.name)
        print("[MicroHapulator::sim]", message, file=sys.stderr)
    if args.haplo_seq:
        index = MicrohapIndex.from_files(args.markers, fasta_path=args.sequences)
        index.validate(symmetric=True)
        with open(args.haplo_seq, "w") as fh:
            for defline, sequence in profile.haploseqs(index):
                print(">", defline, "\n", sequence, sep="", file=fh)
            message = "haplotype sequences written to {:s}".format(fh.name)
            print("[MicroHapulator::sim]", message, file=sys.stderr)
