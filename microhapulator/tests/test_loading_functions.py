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


@pytest.mark.parametrize(
    "tsv,nrows,markerid,offset",
    [
        ("def/acb-dozen-offsets.tsv", 42, "mh09USC-9qB", 193),
        ("def/yellow-offsets.tsv", 21, "mh01KK-205", 207),
    ],
)
def test_load_marker_definitions(tsv, nrows, markerid, offset):
    markers = microhapulator.load_marker_definitions(data_file(tsv))
    assert list(markers.columns) == ["Marker", "Offset"]
    assert markers.shape[0] == nrows
    assert markers.Marker.iloc[10] == markerid
    assert markers.Offset.iloc[10] == offset


def test_load_marker_definitions_missing_column():
    message = r"column\(s\) missing from marker definition file: Marker"
    with pytest.raises(ValueError, match=message):
        microhapulator.load_marker_definitions(data_file("def/orange-offsets-missing.tsv"))


@pytest.mark.parametrize(
    "fasta,nseq,firstseq",
    [
        ("refr/acb-dozen-refr.fasta", 12, "mh07USC-7qB"),
        ("refr/russ4-refr.fasta.gz", 4, "mh01KK-001"),
    ],
)
def test_load_marker_reference_sequences(fasta, nseq, firstseq):
    seqdict = microhapulator.load_marker_reference_sequences(data_file(fasta))
    assert len(seqdict) == nseq
    observed = next(iter(seqdict.keys()))
    print(firstseq, observed)
    assert firstseq == observed


def test_cross_check_equal():
    s1 = set(["mh19KK-056", "mh07NH-12", "mh06USC-6qD"])
    s2 = set(["mh19KK-056", "mh07NH-12", "mh06USC-6qD"])
    microhapulator.cross_check_marker_ids(s1, s2, "marker definitions", "marker frequencies")


def test_cross_check_unequal():
    s1 = set(["mh19KK-056", "mh07NH-12", "mh06USC-6qD"])
    s2 = set(["mh19KK-056", "mh07NH-12"])
    with pytest.raises(ValueError, match=r"marker IDs unique to set1={mh06USC-6qD};"):
        microhapulator.cross_check_marker_ids(s1, s2, "marker definitions", "marker frequencies")
    with pytest.raises(ValueError, match=r"marker IDs unique to set2={mh06USC-6qD};"):
        microhapulator.cross_check_marker_ids(s2, s1, "marker definitions", "marker frequencies")
