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
import re


@dataclass
class FilterStats:
    filtered: int
    retained: int

    @property
    def total(self):
        return self.filtered + self.retained

    @property
    def total_reads(self):
        return f"{self.total:,}"

    @property
    def excluded(self):
        return f"{self.filtered:,}"

    @property
    def kept(self):
        return f"{self.retained:,}"

    @property
    def retention_rate(self):
        rate = self.retained / self.total * 100
        return f"{rate:.1f}%"

    @classmethod
    def from_txt(cls, filepath):
        table = pd.read_csv(filepath, sep="\t")
        assert len(table) == 1
        entry = table.iloc[0]
        stats = cls(entry.ReadsRemoved, entry.ReadsKept)
        return stats


@dataclass
class PairedAmbiguityFilterStats:
    total: int
    failed_r1: int
    failed_r2: int
    failed_both: int
    retained: int

    @property
    def total_reads(self):
        return f"{self.total:,}"

    @property
    def excluded(self):
        failed = self.failed_r1 + self.failed_r2 + self.failed_both
        return f"{failed:,}"

    @property
    def excluded_r1(self):
        return f"{self.failed_r1:,}"

    @property
    def excluded_r2(self):
        return f"{self.failed_r2:,}"

    @property
    def excluded_both(self):
        return f"{self.failed_both:,}"

    @property
    def kept(self):
        return f"{self.retained:,}"

    @property
    def retention_rate(self):
        rate = self.retained / self.total
        return f"{rate * 100:.1f}%"

    @classmethod
    def from_txt(cls, filepath):
        table = pd.read_csv(filepath, sep="\t")
        assert len(table) == 1
        entry = table.iloc[0]
        stats = cls(
            entry.TotalPairs,
            entry.R1OnlyFailed,
            entry.R2OnlyFailed,
            entry.PairsFailed,
            entry.PairsKept,
        )
        return stats


@dataclass
class ReadMergingStats:
    total: int
    combined: int
    innies: int
    outies: int
    uncombined: int

    @property
    def total_reads(self):
        return f"{self.total:,}"

    @property
    def merged_reads(self):
        return f"{self.combined:,}"

    @property
    def merge_rate(self):
        rate = self.combined / self.total
        return f"{rate * 100:.1f}%"

    @classmethod
    def from_log(cls, logpath):
        with open(logpath, "r") as fh:
            data = fh.read()
            tp = re.search(r"Total pairs:\s+(\d+)", data).group(1)
            cp = re.search(r"Combined pairs:\s+(\d+)", data).group(1)
            ip = re.search(r"Innie pairs:\s+(\d+)", data).group(1)
            op = re.search(r"Outie pairs:\s+(\d+)", data).group(1)
            up = re.search(r"Uncombined pairs:\s+(\d+)", data).group(1)
            tp, cp, ip, op, up = map(int, (tp, cp, ip, op, up))
            assert cp + up == tp, logpath
            stats = cls(tp, cp, ip, op, up)
            return stats
