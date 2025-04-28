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

from collections import namedtuple, defaultdict
from dataclasses import dataclass
import pandas as pd
from pathlib import Path
import re
from typing import Dict


class MappingSummary(dict):
    @classmethod
    def from_workdir(cls, samples, workdir="."):
        summary = cls()
        for sample in sorted(samples):
            summary[sample] = MappingStats.from_workdir(sample, workdir=workdir)
        return summary

    @property
    def markers(self):
        for stats in self.values():
            return sorted(stats.repetitive.keys())

    def repetitive_reads_by_marker(self):
        counts = defaultdict(dict)
        for marker in self.markers:
            for sample, stats in self.items():
                counts[marker][sample] = stats.read_count_pair(marker)
        return counts


@dataclass
class MappingStats:
    total: int
    mapped: int
    chisq: float
    marker_counts: Dict[str, int]
    repetitive: Dict[str, int]

    @property
    def total_reads(self):
        return f"{self.total:,}"

    @property
    def mapped_reads(self):
        return f"{self.mapped:,}"

    @property
    def mapping_rate(self):
        rate = self.mapped / self.total * 100
        return f"{rate:.1f}%"

    @property
    def chi_square(self):
        return f"{self.chisq:.3f}"

    def read_count_pair(self, marker):
        if marker not in self.marker_counts or marker not in self.repetitive:
            return KeyError(marker)
        return ReadCountPair(self.marker_counts[marker], self.repetitive[marker])

    @classmethod
    def from_workdir(cls, sample, workdir="."):
        align_dir = Path(workdir) / "analysis" / sample / "alignment"
        stats_file = align_dir / f"{sample}.bam.stats"
        chisq_file = align_dir / f"{sample}-interlocus-balance-chisq.txt"
        rep_file = align_dir / f"{sample}-repetitive-reads.csv"
        counts_file = align_dir / f"{sample}-marker-read-counts.csv"
        total_reads, mapped_reads = MappingStats.parse_read_stats(stats_file)
        chi_square = MappingStats.parse_chi_square_stat(chisq_file)
        map_read_counts = pd.read_csv(counts_file).set_index("Marker")["ReadCount"].to_dict()
        rep_read_counts = pd.read_csv(rep_file).set_index("Marker")["RepetitiveReads"].to_dict()
        return cls(total_reads, mapped_reads, chi_square, map_read_counts, rep_read_counts)

    @staticmethod
    def parse_read_stats(stats_file):
        with open(stats_file, "r") as fh:
            stats = fh.read()
            total_reads = re.search("SN\traw total sequences:\t(\d+)", stats).group(1)
            mapped_reads = re.search("SN\treads mapped:\t(\d+)", stats).group(1)
            return int(total_reads), int(mapped_reads)

    @staticmethod
    def parse_chi_square_stat(chisq_file):
        with open(chisq_file, "r") as fh:
            chisq = fh.read()
            chi_square = re.search("\(chi-square statistic\): (\S+)", chisq).group(1)
            return float(chi_square)


ReadCountPair = namedtuple("CountPair", "mapped,repetitive")
