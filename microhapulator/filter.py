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

from Bio import SeqIO
from dataclasses import dataclass
from microhapulator import open as mhopen


class PairedReadFilter:
    def __init__(self, r1_in, r2_in, output_prefix):
        self.r1_in = mhopen(r1_in, "r")
        self.r2_in = mhopen(r2_in, "r")
        self.outstream = PairedReadOutputStream(output_prefix)
        self.counts = PairedReadFilterCounts()

    def __iter__(self):
        with self.r1_in as r1_in, self.r2_in as r2_in:
            r1 = SeqIO.parse(r1_in, "fastq")
            r2 = SeqIO.parse(r2_in, "fastq")
            for r1, r2 in zip(r1, r2):
                yield r1, r2

    def filter(self):
        for r1, r2 in self:
            keep_r1, keep_r2 = self.keep(r1, r2)
            if keep_r1 and keep_r2:
                SeqIO.write(r1, self.outstream.r1, "fastq")
                SeqIO.write(r2, self.outstream.r2, "fastq")
                self.counts.pairs_passed += 1
            elif keep_r1:
                SeqIO.write(r1, self.outstream.r2_mates, "fastq")
                self.counts.r2_failed += 1
            elif keep_r2:
                SeqIO.write(r2, self.outstream.r1_mates, "fastq")
                self.counts.r1_failed += 1
            else:
                self.counts.pairs_failed += 1


@dataclass
class PairedReadFilterCounts:
    r1_failed: int = 0
    r2_failed: int = 0
    pairs_failed: int = 0
    pairs_passed: int = 0

    @property
    def total_pairs(self):
        return self.total_pairs_filtered + self.pairs_passed

    @property
    def total_pairs_filtered(self):
        return self.r1_failed + self.r2_failed + self.pairs_failed


class PairedReadOutputStream:
    def __init__(self, prefix):
        self.r1 = open(f"{prefix}-ambig-filtered-R1.fastq", "w")
        self.r2 = open(f"{prefix}-ambig-filtered-R2.fastq", "w")
        self.r1_mates = open(f"{prefix}-ambig-R1-mates.fastq", "w")
        self.r2_mates = open(f"{prefix}-ambig-R2-mates.fastq", "w")

    def __del__(self):
        self.r1.close()
        self.r2.close()
        self.r1_mates.close()
        self.r2_mates.close()


class AmbigPairedReadFilter(PairedReadFilter):
    def __init__(self, r1_in, r2_in, out_prefix, threshold=0.2):
        self.threshold = threshold
        super().__init__(r1_in, r2_in, out_prefix)

    def keep(self, r1, r2):
        keep_r1 = not self.is_ambiguous(r1)
        keep_r2 = not self.is_ambiguous(r2)
        return keep_r1, keep_r2

    def is_ambiguous(self, read):
        return read.seq.count("N") / len(read.seq) > self.threshold

    @property
    def summary(self):
        header = "R1OnlyFailed\tR2OnlyFailed\tPairsFailed\tPairsRemoved\tPairsKept\tTotalPairs"
        count_data = f"{self.counts.r1_failed}\t{self.counts.r2_failed}\t{self.counts.pairs_failed}\t{self.counts.total_pairs_filtered}\t{self.counts.pairs_passed}\t{self.counts.total_pairs}"
        return f"{header}\n{count_data}"


class SingleReadFilter:
    def __init__(self, reads_in, reads_out):
        self.reads_in = mhopen(
            reads_in,
            "r",
        )
        self.reads_out = open(reads_out, "w")
        self.num_reads_failed = 0
        self.num_reads_passed = 0

    def __del__(self):
        self.reads_out.close()

    def filter(self):
        with self.reads_in as r_in:
            for read in SeqIO.parse(r_in, "fastq"):
                if self.keep(read):
                    SeqIO.write(read, self.reads_out, "fastq")
                    self.num_reads_passed += 1
                else:
                    self.num_reads_failed += 1


class AmbigSingleReadFilter(SingleReadFilter):
    def __init__(self, reads_in, outfile, threshold=0.2):
        self.threshold = threshold
        super().__init__(reads_in, outfile)

    def keep(self, read):
        is_ambiguous = read.seq.count("N") / len(read.seq) > self.threshold
        return not is_ambiguous

    @property
    def summary(self):
        return f"ReadsRemoved\tReadsKept\n{self.num_reads_failed}\t{self.num_reads_passed}"
