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

from dataclasses import dataclass
import pandas as pd
from typing import Dict


class TypingSummary(dict):
    @classmethod
    def from_workdir(cls, samples, workdir="."):
        summary = cls()
        for sample in sorted(samples):
            summary[sample] = TypingStats.from_workdir(sample, workdir=workdir)
        return summary

    @property
    def has_high_discard_markers(self):
        for stats in self.values():
            if len(stats.high_discard_rates) > 0:
                return True
        return False

    @property
    def high_discard_markers(self):
        tables = list()
        for sample, stats in self.items():
            table = stats.high_discard_rates
            table["Sample"] = sample
            tables.append(table)
        table = pd.concat(tables)
        return table.sort_values(["Marker", "Sample"])

    @property
    def has_high_gap_markers(self):
        for stats in self.values():
            if len(stats.high_gap_rates) > 0:
                return True
        return False

    @property
    def high_gap_markers(self):
        tables = list()
        for sample, stats in self.items():
            table = stats.high_gap_rates
            table["Sample"] = sample
            tables.append(table)
        table = pd.concat(tables)
        return table.sort_values(["Marker", "Sample"])


@dataclass
class TypingStats:
    attempted_events: Dict[str, int]
    successful_events: Dict[str, int]
    typing_rates: Dict[str, float]
    high_discard_rates: pd.DataFrame
    high_gap_rates: pd.DataFrame

    @property
    def attempted(self):
        return f"{sum(self.attempted_events.values()):,.0f}"

    @property
    def successful(self):
        return f"{sum(self.successful_events.values()):,.0f}"

    @property
    def typing_rate(self):
        numerator = sum(self.successful_events.values())
        denominator = sum(self.attempted_events.values())
        rate = numerator / denominator * 100
        return f"{rate:.2f}%"

    def marker_total(self, marker):
        if marker not in self.attempted_events or marker not in self.successful_events:
            raise KeyError(marker)
        total = self.attempted_events[marker] + self.successful_events[marker]
        return f"{total:,.0f}"

    def marker_successful(self, marker):
        if marker not in self.successful_events:
            raise KeyError(marker)
        return f"{self.successful_events[marker]:,.0f}"

    def marker_typing_rate(self, marker):
        if marker not in self.typing_rates:
            raise KeyError(marker)
        return f"{self.typing_rates[marker] * 100:.2f}"

    @classmethod
    def from_workdir(cls, sample, workdir="."):
        attempted, successful, rates = {}, {}, {}
        filename = f"{workdir}/analysis/{sample}/03typing/{sample}-typing-rate.tsv"
        table = pd.read_csv(filename, sep="\t").set_index("Marker")
        for marker, row in table.iterrows():
            attempted[marker] = row.TotalReads
            successful[marker] = row.TypedReads
            rates[marker] = row.TypingRate
        filename = f"{workdir}/analysis/{sample}/03typing/{sample}-discard-rate.tsv"
        high_discard = pd.read_csv(filename, sep="\t")
        filename = f"{workdir}/analysis/{sample}/03typing/{sample}-gapped-rate.tsv"
        high_gap = pd.read_csv(filename, sep="\t")
        return cls(attempted, successful, rates, high_discard, high_gap)
