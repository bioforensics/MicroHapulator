# -------------------------------------------------------------------------------------------------
# Copyright (c) 2021, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import microhapulator
from microhapulator.tests import data_file
import pytest


@pytest.mark.parametrize(
    "tsv,nrows,value",
    [
        ("freq/russ4-freq.tsv", 38, 0.123),
        ("freq/acb-dozen-freq.tsv", 61, 0.182),
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


def test_load_marker_filters():
    markers = ("mh01XYZ-1", "mh02XYZ-2", "mh02XYZ-3")
    configfile = data_file("filters.csv")
    thresholds = microhapulator.load_marker_thresholds(
        markers, configfile=configfile, global_static=10, global_dynamic=0.015
    )
    print(thresholds)
    assert thresholds.shape == (3, 2)
    assert thresholds.loc["mh01XYZ-1", "Static"] == 5
    assert thresholds.loc["mh02XYZ-2", "Static"] == 10
    assert thresholds.loc["mh02XYZ-3", "Static"] == 10
    assert thresholds.loc["mh01XYZ-1", "Dynamic"] == pytest.approx(0.001)
    assert thresholds.loc["mh02XYZ-2", "Dynamic"] == pytest.approx(0.015)
    assert thresholds.loc["mh02XYZ-3", "Dynamic"] == pytest.approx(0.015)


def test_load_marker_filters_default():
    markers = ("mh04ZZZ-1", "mh06ZZZ-2", "mh19ZZZ-3", "mh22ZZZ-4")
    thresholds = microhapulator.load_marker_thresholds(markers)
    assert thresholds.shape == (4, 2)
    for marker, row in thresholds.iterrows():
        assert row.Static == 5
        assert row.Dynamic == pytest.approx(0.02)
