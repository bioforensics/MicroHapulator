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

import microhapulator
from microhapulator.op import sim
import sys


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


def load_inputs(freqfile, markerfile, seqfile, haploseqs=False):
    frequencies = microhapulator.load_marker_frequencies(freqfile)
    if not haploseqs:
        return frequencies, None, None
    markers = microhapulator.load_marker_definitions(markerfile)
    sequences = microhapulator.load_marker_reference_sequences(seqfile)
    microhapulator.cross_check_marker_ids(
        frequencies.Marker, markers.Marker, "marker frequencies", "marker definitions"
    )
    microhapulator.cross_check_marker_ids(
        frequencies.Marker, sequences.keys(), "marker frequencies", "marker reference sequences"
    )
    return frequencies, markers, sequences


def main(args):
    frequencies, markers, sequences = load_inputs(
        args.freq, args.markers, args.sequences, haploseqs=args.haplo_seq
    )
    profile = sim(frequencies, seed=args.seed)
    with microhapulator.open(args.out, "w") as fh:
        profile.dump(fh)
        message = "profile JSON written to {:s}".format(fh.name)
        print("[MicroHapulator::sim]", message, file=sys.stderr)
    if args.haplo_seq:
        with microhapulator.open(args.haplo_seq, "w") as fh:
            for defline, sequence in profile.haploseqs(markers, sequences):
                print(">", defline, "\n", sequence, sep="", file=fh)
            message = "haplotype sequences written to {:s}".format(fh.name)
            print("[MicroHapulator::sim]", message, file=sys.stderr)
