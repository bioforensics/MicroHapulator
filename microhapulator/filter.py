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
from microhapulator import open as mhopen
import gzip


class PairedReadFilter:
    def __init__(self, r1_in, r2_in):
        self.r1_in = mhopen(r1_in, "r")
        self.r2_in = mhopen(r2_in, "r")
        self.r1_out = mhopen(self.r1_out, "w")
        self.r2_out = mhopen(self.r2_out, "w")
        self.r1_mates_out = mhopen(self.r1_mates_out, "w")
        self.r2_mates_out = mhopen(self.r2_mates_out, "w")
        self.counts_out = mhopen(self.counts_out, "w")
        self.num_r1_failed = 0
        self.num_r2_failed = 0
        self.num_both_failed = 0

    def __iter__(self):
        with self.r1_in as r1_in, self.r2_in as r2_in:
            r1 = SeqIO.parse(r1_in, "fastq")
            r2 = SeqIO.parse(r2_in, "fastq")
            for r1, r2 in zip(r1, r2):
                yield r1, r2

    def filter(self):
        with self.r1_out as r1_out, self.r2_out as r2_out, self.r1_mates_out as r1_mates, self.r2_mates_out as r2_mates:
            for r1, r2 in self:
                keep_r1, keep_r2 = self.keep(r1, r2)
                if keep_r1 and keep_r2:
                    SeqIO.write(r1, r1_out, "fastq")
                    SeqIO.write(r2, r2_out, "fastq")
                elif keep_r1:
                    SeqIO.write(r1, r2_mates, "fastq")
                    self.num_r2_failed += 1
                elif keep_r2:
                    SeqIO.write(r2, r1_mates, "fastq")
                    self.num_r1_failed += 1
                else:
                    self.num_both_failed += 1
        self.write_counts_output()


class AmbigPairedReadFilter(PairedReadFilter):
    def __init__(self, r1_in, r2_in, out_prefix, threshold=0.2):
        self.r1_out = f"{out_prefix}-ambig-filtered-R1.fastq"
        self.r2_out = f"{out_prefix}-ambig-filtered-R2.fastq"
        self.r1_mates_out = f"{out_prefix}-ambig-R1-mates.fastq"
        self.r2_mates_out = f"{out_prefix}-ambig-R2-mates.fastq"
        self.counts_out = f"{out_prefix}-ambig-read-counts.txt"
        self.threshold = threshold
        super().__init__(r1_in, r2_in)

    def keep(self, r1, r2):
        keep_r1 = not self.is_ambiguous(r1)
        keep_r2 = not self.is_ambiguous(r2)
        return keep_r1, keep_r2

    def is_ambiguous(self, read):
        return read.count("N") / len(read) > self.threshold

    def write_counts_output(self):
        with self.counts_out as counts_out:
            header = f"R1Only\tR2Only\tR1andR2\tPairsRemoved"
            total_pairs_filtered = self.num_r1_failed + self.num_r2_failed + self.num_both_failed
            counts_str = f"{self.num_r1_failed}\t{self.num_r2_failed}\t{self.num_both_failed}\t{total_pairs_filtered}"
            counts_out.write(f"{header}\n{counts_str}\n")


class SingleReadFilter:
    def __init__(self, reads_in):
        self.reads_in = mhopen(reads_in, "r")
        self.reads_out = mhopen(self.reads_out, "w")
        self.counts_out = mhopen(self.counts_out, "w")
        self.num_reads_failed = 0

    def filter(self):
        with self.reads_in as r_in, self.reads_out as r_out:
            for read in SeqIO.parse(r_in, "fastq"):
                if self.keep(read):
                    SeqIO.write(read, r_out, "fastq")
                else:
                    self.num_reads_failed += 1
        self.write_counts_output()


class AmbigSingleReadFilter(SingleReadFilter):
    def __init__(self, reads_in, output_prefix, threshold=0.2):
        self.counts_out = f"{output_prefix}-ambig-read-counts.txt"
        self.reads_out = f"{output_prefix}-ambig-filtered.fastq"
        self.threshold = threshold
        super().__init__(reads_in)

    def keep(self, read):
        is_ambiguous = read.count("N") / len(read) > self.threshold
        return not is_ambiguous

    def write_counts_output(self):
        with self.counts_out as counts_out:
            counts_out.write(f"ReadsRemoved\n{self.num_reads_failed}\n")
