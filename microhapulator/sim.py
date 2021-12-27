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
from microhapulator.profile import SimulatedProfile
import numpy.random


def sim(frequencies, seed=None):
    """Simulate a diploid genotype from the specified microhaplotype frequencies."""
    profile = SimulatedProfile(ploidy=2)
    if seed is None:
        seed = numpy.random.randint(2 ** 32 - 1)
    profile.data["metadata"] = {
        "HaploSeed": seed,
    }
    numpy.random.seed(seed)
    markers = sorted(frequencies.Marker.unique())
    for haploindex in range(2):
        for marker in markers:
            haplofreqs = frequencies[frequencies.Marker == marker]
            haplotypes = list(haplofreqs.Haplotype)
            freqs = list(haplofreqs.Frequency)
            freqs = [x / sum(freqs) for x in freqs]
            sampled_haplotype = numpy.random.choice(haplotypes, p=freqs)
            profile.add(haploindex, marker, sampled_haplotype)
    message = f"simulated microhaplotype variation at {len(markers)} markers"
    microhapulator.plog("[MicroHapulator::sim]", message)
    return profile


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
        microhapulator.plog("[MicroHapulator::sim]", message)
    if args.haplo_seq:
        with microhapulator.open(args.haplo_seq, "w") as fh:
            for defline, sequence in profile.haploseqs(markers, sequences):
                print(">", defline, "\n", sequence, sep="", file=fh)
            message = "haplotype sequences written to {:s}".format(fh.name)
            microhapulator.plog("[MicroHapulator::sim]", message)
