import microhapulator.api as mhapi
from microhapulator.tests import data_file
import pandas as pd
import pytest


@pytest.mark.parametrize("ambig_thresh,expected", [(0.2, True), (0.5, False)])
def test_is_ambiguous(ambig_thresh, expected):
    sequence = "ANTACNTANGACNCGATACAN"
    is_ambiguous = mhapi.is_ambiguous(sequence, ambiguous_thresh=ambig_thresh)
    assert is_ambiguous == expected


@pytest.mark.parametrize(
    "fastq", [data_file("ambiguous-single-end.fastq"), data_file("ambiguous-single-end.fastq.gz")]
)
def test_filter_ambiguous_single_reads(tmp_path, fastq):
    mhapi.filter_ambiguous_single_reads(fastq, tmp_path, sample="test-sample")
    filtered_reads_file = tmp_path / "test-sample-ambiguous-filtered.fastq"
    counts = tmp_path / "test-sample-ambiguous-read-counts.txt"
    assert filtered_reads_file.is_file()
    assert counts.is_file()
    counts_df = pd.read_csv(counts, sep="\t")
    expected_reads_removed = 2
    assert counts_df["ReadsRemoved"][0] == expected_reads_removed
    filtered_reads = mhapi.parse_reads(str(filtered_reads_file))
    assert len(filtered_reads) == 8


@pytest.mark.parametrize(
    "r1,r2", [(data_file("ambiguous-r1.fastq"), data_file("ambiguous-r2.fastq"))]
)
def test_filter_ambiguous_paried_reads(tmp_path, r1, r2):
    mhapi.filter_ambiguous_paired_reads(r1, r2, tmp_path, sample="test-sample")
    filtered_r1_file = tmp_path / "test-sample-r1-ambiguous-filtered.fastq"
    filtered_r2_file = tmp_path / "test-sample-r2-ambiguous-filtered.fastq"
    r1_mates_file = tmp_path / "test-sample-ambiguous-r1-mates.fastq"
    r2_mates_file = tmp_path / "test-sample-ambiguous-r2-mates.fastq"
    counts = tmp_path / "test-sample-ambiguous-read-counts.txt"
    assert filtered_r1_file.is_file()
    assert filtered_r2_file.is_file()
    assert r1_mates_file.is_file()
    assert r2_mates_file.is_file()
    assert counts.is_file()
    assert len(mhapi.parse_reads(str(filtered_r1_file))) == 7
    assert len(mhapi.parse_reads(str(filtered_r2_file))) == 7
    assert len(mhapi.parse_reads(str(r1_mates_file))) == 1
    assert len(mhapi.parse_reads(str(r2_mates_file))) == 1
    counts_df = pd.read_csv(counts, sep="\t")
    expected_counts_df = pd.read_csv(data_file("ambiguous-read-counts-paired.txt"), sep="\t")
    assert counts_df.equals(expected_counts_df)
