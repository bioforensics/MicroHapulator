from Bio import SeqIO
from microhapulator.filter import AmbigPairedReadFilter, AmbigSingleReadFilter
from microhapulator.tests import data_file
import pandas as pd
import pytest


def test_filter_ambiguous_single_reads(tmp_path):
    reads = data_file("ambiguous-single-end.fastq")
    ambig_filter = AmbigSingleReadFilter(reads, tmp_path / "test")
    ambig_filter.filter()
    ambig_filter.write_counts_output()
    filtered_reads_file = tmp_path / "test-ambig-filtered.fastq"
    counts = tmp_path / "test-ambig-read-counts.txt"
    assert filtered_reads_file.is_file()
    assert counts.is_file()
    counts_df = pd.read_csv(counts, sep="\t")
    expected_reads_removed = 2
    assert counts_df["ReadsRemoved"][0] == expected_reads_removed
    assert len(list(SeqIO.parse(filtered_reads_file, "fastq"))) == 8


def test_filter_ambiguous_paired_reads(tmp_path):
    r1 = data_file("ambiguous-r1.fastq")
    r2 = data_file("ambiguous-r2.fastq")
    ambig_filter = AmbigPairedReadFilter(r1, r2, tmp_path / "test")
    ambig_filter.filter()
    ambig_filter.write_counts_output()
    filtered_r1_file = tmp_path / "test-ambig-filtered-R1.fastq"
    filtered_r2_file = tmp_path / "test-ambig-filtered-R2.fastq"
    r1_mates_file = tmp_path / "test-ambig-R1-mates.fastq"
    r2_mates_file = tmp_path / "test-ambig-R2-mates.fastq"
    counts = tmp_path / "test-ambig-read-counts.txt"
    assert filtered_r1_file.is_file()
    assert filtered_r2_file.is_file()
    assert r1_mates_file.is_file()
    assert r2_mates_file.is_file()
    assert counts.is_file()
    assert len(list(SeqIO.parse(filtered_r1_file, "fastq"))) == 6
    assert len(list(SeqIO.parse(filtered_r2_file, "fastq"))) == 6
    assert len(list(SeqIO.parse(r1_mates_file, "fastq"))) == 2
    assert len(list(SeqIO.parse(r2_mates_file, "fastq"))) == 1
    counts_df = pd.read_csv(counts, sep="\t")
    expected_counts_df = pd.read_csv(data_file("ambiguous-read-counts-paired.txt"), sep="\t")
    assert counts_df.equals(expected_counts_df)
