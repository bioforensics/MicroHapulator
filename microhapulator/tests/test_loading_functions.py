# -----------------------------------------------------------------------------
# Copyright (c) 2021, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
from microhapulator.tests import data_file
import pytest


@pytest.mark.parametrize(
    "tsv,nrows,value",
    [
        ("freq/russ4-freq.tsv", 38, 0.123),
        ("acb-dozen-freq.tsv", 61, 0.182),
    ],
)
def test_load_marker_frequencies(tsv, nrows, value):
    freqs = microhapulator.load_marker_frequencies(data_file(tsv))
    assert list(freqs.columns) == ["Marker", "Haplotype", "Frequency"]
    assert freqs.shape[0] == nrows
    assert freqs.Frequency.iloc[23] == pytest.approx(value)


def test_load_marker_frequencies_extra_column_ok():
    freqs = microhapulator.load_marker_frequencies(data_file("freq/korea-5loc-freq-extracol.tsv"))
    assert list(freqs.columns) == ["Marker", "Haplotype", "Frequency", "Foo"]
    assert freqs.Haplotype.iloc[1] == "C,A,G,G"


@pytest.mark.parametrize(
    "tsv",
    [
        "freq/korea-5loc-freq-missingcol.tsv",
        "freq/korea-5loc-freq-badcol.tsv",
    ],
)
def test_load_marker_frequencies_missing_or_bad_columns(tsv):
    message = r"column\(s\) missing from marker frequency file: Haplotype"
    with pytest.raises(ValueError, match=message):
        microhapulator.load_marker_frequencies(data_file(tsv))
