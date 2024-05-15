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

from .qcstats import FilterStats, PairedAmbiguityFilterStats, ReadMergingStats
from dataclasses import dataclass


@dataclass
class SingleEndReadQCSummary:
    ambig: FilterStats
    length: FilterStats

    @property
    def total(self):
        return self.ambig.filtered + self.ambig.retained

    @property
    def total_reads(self):
        return f"{self.total:,}"

    @property
    def filtered_ambig(self):
        return SingleEndReadQCSummary.render(self.ambig.filtered, self.total)

    @property
    def filtered_length(self):
        return SingleEndReadQCSummary.render(self.length.filtered, self.total)

    @property
    def retention(self):
        return SingleEndReadQCSummary.render(self.length.retained, self.total)

    @staticmethod
    def render(numerator, denominator):
        rate = numerator / denominator * 100
        return f"{numerator:,} ({rate:.1f}%)"

    @classmethod
    def from_workdir(cls, wdpath, sample):
        ambig_path = f"{wdpath}/analysis/{sample}/{sample}-ambig-read-counts.txt"
        length_path = f"{wdpath}/analysis/{sample}/{sample}-length-filtered-read-counts.txt"
        ambig_stats = FilterStats.from_txt(ambig_path)
        length_stats = FilterStats.from_txt(length_path)
        return cls(ambig_stats, length_stats)

    @staticmethod
    def collect(samples, workdir="."):
        sampleqc = dict()
        for sample in sorted(samples):
            sampleqc[sample] = SingleEndReadQCSummary.from_workdir(workdir, sample)
        return sampleqc


@dataclass
class PairedReadQCSummary:
    ambig: PairedAmbiguityFilterStats
    merge: ReadMergingStats
    length: FilterStats

    @classmethod
    def from_workdir(cls, wdpath, sample):
        ambig_path = f"{wdpath}/analysis/{sample}/{sample}-ambig-read-counts.txt"
        merge_path = f"{wdpath}/analysis/{sample}/flash.log"
        length_path = f"{wdpath}/analysis/{sample}/{sample}-length-filtered-read-counts.txt"
        ambig_stats = PairedAmbiguityFilterStats.from_txt(ambig_path)
        merge_stats = ReadMergingStats.from_log(merge_path)
        length_stats = FilterStats.from_txt(length_path)
        return cls(ambig_stats, merge_stats, length_stats)

    @staticmethod
    def collect(samples, workdir="."):
        sampleqc = dict()
        for sample in sorted(samples):
            sampleqc[sample] = PairedReadQCSummary.from_workdir(workdir, sample)
        return sampleqc
