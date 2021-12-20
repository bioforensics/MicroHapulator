#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from Bio import SeqIO
from collections import defaultdict
from happer.mutate import mutate
import microhapulator
from microhapulator.profile import SimulatedProfile
import numpy.random
import pandas as pd


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
    frequencies = pd.read_csv(freqfile, sep="\t")
    columns = list(frequencies.columns)
    assert "Marker" in columns
    assert "Haplotype" in columns
    assert "Frequency" in columns
    if not haploseqs:
        return frequencies, None, None
    markers = pd.read_csv(markerfile, sep="\t")
    columns = list(markers.columns)
    assert "Marker" in columns
    assert "Offset" in columns
    sequences = SeqIO.to_dict(SeqIO.parse(seqfile, "fasta"))
    sequences = {seqid: record.seq for seqid, record in sequences.items()}
    if set(frequencies.Marker) != set(markers.Marker):
        frequniq = set(frequencies.Marker) - set(markers.Marker)
        markeruniq = set(markers.Marker) - set(frequencies.Marker)
        message = "discrepancy between marker definitions and population frequencies:"
        if frequniq:
            markerids = ", ".join(frequniq)
            message += f" markers with frequency data but no definition={{{markerids}}};"
        if markeruniq:
            markerids = ", ".join(markeruniq)
            message += f" markers with defined offsets but no frequency data={{{markerids}}};"
        raise ValueError(message)
    if set(frequencies.Marker) != set(sequences.keys()):
        frequniq = set(frequencies.Marker) - set(sequences.keys())
        sequniq = set(sequences.keys()) - set(frequencies.Marker)
        message = "discrepancy between marker definitions and population frequencies:"
        if frequniq:
            markerids = ", ".join(frequniq)
            message += f" markers with frequency data but no reference sequence={{{markerids}}};"
        if sequniq:
            markerids = ", ".join(sequniq)
            message += f" markers with a reference sequence but no frequency data={{{markerids}}};"
        raise ValueError(message)
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
