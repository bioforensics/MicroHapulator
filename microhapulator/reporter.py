# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from .mapstats import MappingSummary
from .qcsummary import PairedReadQCSummary, SingleEndReadQCSummary
import json
import pandas as pd
from pathlib import Path


class Reporter:
    def __init__(self, samples, workdir=".", reads_are_paired=True):
        if reads_are_paired:
            self.read_length_table = read_length_table_paired_end(samples, workdir=workdir)
            self.qc_summary = PairedReadQCSummary.collect(samples, workdir=workdir)
        else:
            self.read_length_table = read_length_table_single_end(samples, workdir=workdir)
            self.qc_summary = SingleEndReadQCSummary.collect(samples, workdir=workdir)
        self.mapping_summary = MappingSummary.from_workdir(samples)
        self.per_marker_typing_rates = per_marker_typing_rate(samples, workdir=workdir)
        self.per_marker_mapping_rates = per_marker_mapping_rate(samples, workdir=workdir)
        self.thresholds = None

    @property
    def marker_names(self):
        for sample_rates in self.per_marker_mapping_rates.values():
            return sample_rates.index


def read_length_table_paired_end(samples, workdir="."):
    read_length_data = list()
    for sample in samples:
        with open(f"{workdir}/analysis/{sample}/{sample}-r1-read-lengths.json", "r") as fh:
            r1lengths = json.load(fh)
            r1lengths = list(set(r1lengths))
        with open(f"{workdir}/analysis/{sample}/{sample}-r2-read-lengths.json", "r") as fh:
            r2lengths = json.load(fh)
            r2lengths = list(set(r2lengths))
        if len(r1lengths) != 1 or len(r2lengths) != 1:
            return None
        read_length_data.append((sample, r1lengths[0], r2lengths[0]))
    else:
        return pd.DataFrame(read_length_data, columns=("Sample", "LengthR1", "LengthR2"))


def read_length_table_single_end(samples, workdir="."):
    read_length_data = list()
    for sample in samples:
        with open(f"{workdir}/analysis/{sample}/{sample}-read-lengths.json", "r") as fh:
            lengths = json.load(fh)
            lengths = list(set(lengths))
        if len(lengths) != 1:
            return None
        read_length_data.append((sample, lengths[0]))
    return pd.DataFrame(read_length_data, columns=("Sample", "Length"))


def per_marker_typing_rate(samples, workdir="."):
    sample_rates = dict()
    for sample in samples:
        filename = f"{workdir}/analysis/{sample}/{sample}-typing-rate.tsv"
        sample_df = pd.read_csv(filename, sep="\t").set_index("Marker")
        sample_rates[sample] = sample_df
    return sample_rates


def per_marker_mapping_rate(samples, workdir="."):
    sample_rates = dict()
    for sample in samples:
        total_reads_filename = f"{workdir}/analysis/{sample}/{sample}-marker-read-counts.csv"
        total_df = pd.read_csv(total_reads_filename).set_index("Marker")
        expected_count = total_df["ReadCount"].sum() / len(total_df)
        total_df["ExpectedObservedRatio"] = round(total_df["ReadCount"] / expected_count, 2)
        repetitive_filename = Path(f"{workdir}/analysis/{sample}/{sample}-repetitive-reads.csv")
        if repetitive_filename.stat().st_size > 0:
            repetitive_df = pd.read_csv(repetitive_filename).set_index("Marker")
            sample_df = pd.concat([total_df, repetitive_df], axis=1)
            sample_df["RepetitiveRate"] = sample_df["RepetitiveReads"] / sample_df["ReadCount"]
        else:
            sample_df = total_df
            sample_df["RepetitiveReads"] = None
            sample_df["RepetitiveRate"] = None
        sample_rates[sample] = sample_df
    # marker_names = sample_df.index #FIXME
    return sample_rates
