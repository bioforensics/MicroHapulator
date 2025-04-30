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

from .details import MarkerDetails
from .mapstats import MappingSummary
from .qcsummary import PairedReadQCSummary, SingleEndReadQCSummary
from .typestats import TypingSummary
from datetime import datetime
from importlib.resources import files
from jinja2 import FileSystemLoader, Environment, Template
import microhapulator
from microhapulator.marker import MicrohapIndex
import json
import pandas as pd
from pathlib import Path


class OverviewReporter:
    def __init__(self, samples, thresholds, workdir=".", reads_are_paired=True):
        self.samples = sorted(samples)
        self.thresholds = thresholds
        self.reads_are_paired = reads_are_paired
        if reads_are_paired:
            self.read_length_table = read_length_table_paired_end(samples, workdir=workdir)
            self.qc_summary = PairedReadQCSummary.collect(samples, workdir=workdir)
        else:
            self.read_length_table = read_length_table_single_end(samples, workdir=workdir)
            self.qc_summary = SingleEndReadQCSummary.collect(samples, workdir=workdir)
        self.mapping_summary = MappingSummary.from_workdir(samples)
        self.typing_summary = TypingSummary.from_workdir(samples)
        self.per_marker_mapping_rates = per_marker_mapping_rate(samples, workdir=workdir)

    @property
    def marker_names(self):
        for sample_rates in self.per_marker_mapping_rates.values():
            return sample_rates.index

    def render(self):
        template_loader = FileSystemLoader(files("microhapulator") / "data")
        env = Environment(loader=template_loader)
        if self.reads_are_paired:
            template_file = "paired.html"
        else:
            template_file = "single.html"
        template = env.get_template(template_file)
        output = template.render(
            date=datetime.now().replace(microsecond=0).isoformat(),
            mhpl8rversion=microhapulator.__version__,
            samples=self.samples,
            thresholds=self.thresholds,
            read_length_table=self.read_length_table,
            typing_summary=self.typing_summary,
            mapping_rates=self.per_marker_mapping_rates,
            markernames=self.marker_names,
            qc=self.qc_summary,
            mapping_summary=self.mapping_summary,
            repetitive_reads_by_marker=self.mapping_summary.repetitive_reads_by_marker(),
            reads_are_paired=self.reads_are_paired,
            ambiguous_read_threshold=self.thresholds.ambiguous,
            read_length_threshold=self.thresholds.min_read_length,
            discard_alert_threshold=self.thresholds.discard_alert,
            gap_alert_threshold=self.thresholds.gap_alert,
        )
        return output


class DetailReporter:
    def __init__(self, samples, workdir="."):
        self.typing_summary = TypingSummary.from_workdir(samples)
        self.per_marker_mapping_rates = per_marker_mapping_rate(samples, workdir=workdir)
        index = MicrohapIndex.from_files(
            Path(workdir) / "marker-definitions.tsv",
            Path(workdir) / "marker-refr.fasta",
        )
        self.marker_details = list(MarkerDetails.from_index(index))

    @property
    def marker_names(self):
        for sample_rates in self.per_marker_mapping_rates.values():
            return sample_rates.index

    def render(self):
        templatefile = files("microhapulator") / "data" / "marker_details_template.html"
        with open(templatefile, "r") as fh:
            template = Template(fh.read())
            output = template.render(
                date=datetime.now().replace(microsecond=0).isoformat(),
                mhpl8rversion=microhapulator.__version__,
                mapping_rates=self.per_marker_mapping_rates,
                typing_summary=self.typing_summary,
                markernames=sorted(self.marker_names),
                marker_details=self.marker_details,
                isna=pd.isna,
            )
            return output


def read_length_table_paired_end(samples, workdir="."):
    read_length_data = list()
    for sample in samples:
        prefix = f"{workdir}/analysis/{sample}/01preprocessing"
        with open(f"{prefix}/{sample}-r1-read-lengths.json", "r") as fh:
            r1lengths = json.load(fh)
            r1lengths = list(set(r1lengths))
        with open(f"{prefix}/{sample}-r2-read-lengths.json", "r") as fh:
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
        path = f"{workdir}/analysis/{sample}/01preprocessing/{sample}-read-lengths.json"
        with open(path, "r") as fh:
            lengths = json.load(fh)
            lengths = list(set(lengths))
        if len(lengths) != 1:
            return None
        read_length_data.append((sample, lengths[0]))
    return pd.DataFrame(read_length_data, columns=("Sample", "Length"))


def per_marker_mapping_rate(samples, workdir="."):
    sample_rates = dict()
    for sample in samples:
        total_reads_filename = (
            f"{workdir}/analysis/{sample}/02alignment/{sample}-marker-read-counts.csv"
        )
        total_df = pd.read_csv(total_reads_filename).set_index("Marker")
        expected_count = total_df["ReadCount"].sum() / len(total_df)
        total_df["ExpectedObservedRatio"] = total_df["ReadCount"] / expected_count
        repetitive_filename = Path(
            f"{workdir}/analysis/{sample}/02alignment/{sample}-repetitive-reads.csv"
        )
        if repetitive_filename.stat().st_size > 0:
            repetitive_df = pd.read_csv(repetitive_filename).set_index("Marker")
            sample_df = pd.concat([total_df, repetitive_df], axis=1)
            sample_df["RepetitiveRate"] = sample_df["RepetitiveReads"] / sample_df["ReadCount"]
        else:
            sample_df = total_df
            sample_df["RepetitiveReads"] = None
            sample_df["RepetitiveRate"] = None
        sample_rates[sample] = sample_df
    return sample_rates
