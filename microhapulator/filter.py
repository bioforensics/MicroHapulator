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
    def __init__(self, r1_in, r2_in, output_prefix):
        self.r1_in = mhopen(r1_in, "r")
        self.r2_in = mhopen(r2_in, "r")
        self.get_output_files(output_prefix)
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
        with self.outfiles["r1"] as r1_out, self.outfiles["r2"] as r2_out, self.outfiles[
            "r1_mates"
        ] as r1_mates, self.outfiles["r2_mates"] as r2_mates:
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

    def write_counts_output(self):
        with self.outfiles["counts"] as counts_out:
            header = f"R1Only\tR2Only\tR1andR2\tPairsRemoved"
            total_pairs_filtered = self.num_r1_failed + self.num_r2_failed + self.num_both_failed
            counts_str = f"{self.num_r1_failed}\t{self.num_r2_failed}\t{self.num_both_failed}\t{total_pairs_filtered}"
            counts_out.write(f"{header}\n{counts_str}\n")

    def get_output_files(self, out_prefix):
        self.outfiles = dict()
        self.outfiles["r1"] = mhopen(f"{out_prefix}-ambig-filtered-R1.fastq", "w")
        self.outfiles["r2"] = mhopen(f"{out_prefix}-ambig-filtered-R2.fastq", "w")
        self.outfiles["r1_mates"] = mhopen(f"{out_prefix}-ambig-R1-mates.fastq", "w")
        self.outfiles["r2_mates"] = mhopen(f"{out_prefix}-ambig-R2-mates.fastq", "w")
        self.outfiles["counts"] = mhopen(f"{out_prefix}-ambig-read-counts.txt", "w")


class SingleReadFilter:
    def __init__(self, reads_in, output_prefix):
        self.reads_in = mhopen(reads_in, "r")
        self.get_output_files(output_prefix)
        self.num_reads_failed = 0

    def filter(self):
        with self.reads_in as r_in, self.outfiles["reads"] as r_out:
            for read in SeqIO.parse(r_in, "fastq"):
                if self.keep(read):
                    SeqIO.write(read, r_out, "fastq")
                else:
                    self.num_reads_failed += 1


class AmbigSingleReadFilter(SingleReadFilter):
    def __init__(self, reads_in, out_prefix, threshold=0.2):
        self.threshold = threshold
        super().__init__(reads_in, out_prefix)

    def keep(self, read):
        is_ambiguous = read.seq.count("N") / len(read.seq) > self.threshold
        return not is_ambiguous

    def write_counts_output(self):
        with self.outfiles["counts"] as counts_out:
            counts_out.write(f"ReadsRemoved\n{self.num_reads_failed}\n")

    def get_output_files(self, out_prefix):
        self.outfiles = dict()
        self.outfiles["reads"] = mhopen(f"{out_prefix}-ambig-filtered.fastq", "w")
        self.outfiles["counts"] = mhopen(f"{out_prefix}-ambig-read-counts.txt", "w")
