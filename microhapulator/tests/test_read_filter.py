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
from microhapulator.pipe.filter import (
    AmbigPairedReadFilter,
    AmbigSingleReadFilter,
    LengthSingleReadFilter,
)
from microhapulator.tests import data_file
import pandas as pd
import pytest


@pytest.fixture(scope="session")
def ambig_seqs_dir_single(tmp_path_factory):
    tmpdir = tmp_path_factory.mktemp("WD")
    reads = data_file("ambiguous-single-end.fastq")
    ambig_filter = AmbigSingleReadFilter(reads, tmpdir / "test-ambig-filtered.fastq")
    ambig_filter.filter()
    with open(tmpdir / "test-ambig-read-counts.txt", "w") as fh:
        print(ambig_filter.summary, file=fh)
    return tmpdir


@pytest.fixture(scope="session")
def ambig_seqs_dir_paired(tmp_path_factory):
    tmpdir = tmp_path_factory.mktemp("WD")
    r1 = data_file("ambiguous-r1.fastq")
    r2 = data_file("ambiguous-r2.fastq")
    ambig_filter = AmbigPairedReadFilter(r1, r2, tmpdir / "test")
    ambig_filter.filter()
    with open(tmpdir / "test-ambig-read-counts.txt", "w") as fh:
        print(ambig_filter.summary, file=fh)
    return tmpdir


@pytest.fixture(scope="session")
def length_filter_dir(tmp_path_factory):
    tmpdir = tmp_path_factory.mktemp("WD")
    reads = data_file("too-short.fastq")
    length_filter = LengthSingleReadFilter(reads, tmpdir / "test-length-filtered.fastq")
    length_filter.filter()
    with open(tmpdir / "test-length-filtered-read-counts.txt", "w") as fh:
        print(length_filter.summary, file=fh)
    return tmpdir


def test_filter_ambiguous_single_fastq(ambig_seqs_dir_single):
    filtered_reads_file = ambig_seqs_dir_single / "test-ambig-filtered.fastq"
    assert filtered_reads_file.is_file()
    assert len(list(SeqIO.parse(filtered_reads_file, "fastq"))) == 8


def test_filter_ambiguous_single_counts(ambig_seqs_dir_single):
    counts_file = ambig_seqs_dir_single / "test-ambig-read-counts.txt"
    assert counts_file.is_file()
    counts_df = pd.read_csv(counts_file, sep="\t")
    expected_reads_removed = 2
    expected_reads_kept = 8
    assert counts_df["ReadsRemoved"][0] == expected_reads_removed
    assert counts_df["ReadsKept"][0] == expected_reads_kept


def test_filter_ambiguous_paired_fastq(ambig_seqs_dir_paired):
    filtered_r1_file = ambig_seqs_dir_paired / "test-ambig-filtered-R1.fastq"
    filtered_r2_file = ambig_seqs_dir_paired / "test-ambig-filtered-R2.fastq"
    assert filtered_r1_file.is_file()
    assert filtered_r2_file.is_file()
    assert len(list(SeqIO.parse(filtered_r1_file, "fastq"))) == 6
    assert len(list(SeqIO.parse(filtered_r2_file, "fastq"))) == 6


def test_filter_ambiguous_paried_mates(ambig_seqs_dir_paired):
    r1_mates_file = ambig_seqs_dir_paired / "test-ambig-R1-mates.fastq"
    r2_mates_file = ambig_seqs_dir_paired / "test-ambig-R2-mates.fastq"
    assert r1_mates_file.is_file()
    assert r2_mates_file.is_file()
    assert len(list(SeqIO.parse(r1_mates_file, "fastq"))) == 2
    assert len(list(SeqIO.parse(r2_mates_file, "fastq"))) == 1


def test_filter_ambiguous_paired_counts(ambig_seqs_dir_paired):
    counts_file = ambig_seqs_dir_paired / "test-ambig-read-counts.txt"
    assert counts_file.is_file()
    counts_df = pd.read_csv(counts_file, sep="\t")
    expected_counts_df = pd.read_csv(data_file("ambiguous-read-counts-paired.txt"), sep="\t")
    assert counts_df.equals(expected_counts_df)


def test_length_filter_fastq(length_filter_dir):
    filtered_reads_file = length_filter_dir / "test-length-filtered.fastq"
    assert filtered_reads_file.is_file()
    assert len(list(SeqIO.parse(filtered_reads_file, "fastq"))) == 8


def test_length_filter_counts(length_filter_dir):
    counts_file = length_filter_dir / "test-length-filtered-read-counts.txt"
    assert counts_file.is_file()
    counts_df = pd.read_csv(counts_file, sep="\t")
    expected_reads_removed = 2
    expected_reads_kept = 8
    assert counts_df["ReadsRemoved"][0] == expected_reads_removed
    assert counts_df["ReadsKept"][0] == expected_reads_kept
