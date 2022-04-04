# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
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
import microhapulator.api as mhapi
from microhapulator.profile import TypingResult
from microhapulator.tests import data_file
import pytest
from shutil import copyfile


def test_type_simple():
    bam = data_file("pashtun-sim/aligned-reads.bam")
    tsv = data_file("pashtun-sim/tiny-panel.tsv")
    result = mhapi.type(bam, tsv)
    assert result.haplotypes("mh13KK-218") == set()
    assert result.data["markers"]["mh13KK-218"]["typing_result"] == {
        "C,T,C,G": 1,
        "C,T,T,T": 1,
        "G,T,C,T": 1,
        "T,A,C,T": 1,
        "T,G,C,T": 3,
        "T,G,T,T": 2,
        "T,T,A,T": 5,
        "T,T,C,A": 1,
        "T,T,C,C": 2,
        "T,T,C,G": 2,
        "T,T,C,T": 1178,
        "T,T,G,T": 2,
        "T,T,T,A": 2,
        "T,T,T,G": 6,
        "T,T,T,T": 1170,
    }
    assert result.haplotypes("mh21KK-320") == set()
    assert result.data["markers"]["mh21KK-320"]["typing_result"] == {
        "G,A,A,A": 1,
        "G,A,C,A": 3,
        "G,A,G,A": 3,
        "G,A,T,A": 1075,
        "G,A,T,C": 1,
        "G,A,T,G": 1,
        "G,A,T,T": 2,
        "G,C,C,A": 1,
        "G,C,T,A": 4,
        "G,G,A,A": 2,
        "G,G,A,T": 1,
        "G,G,C,A": 1075,
        "G,G,C,C": 3,
        "G,G,C,G": 12,
        "G,G,C,T": 5,
        "G,G,T,A": 4,
        "G,T,C,A": 1,
        "T,G,C,A": 1,
    }


def test_type_missing_bam_index(tmp_path):
    bam = data_file("bam/three-contrib-log-link.bam")
    tsv = data_file("def/default-panel-offsets.tsv")
    tmp_bam = str(tmp_path / "reads.bam")
    tmp_tsv = str(tmp_path / "offsets.tsv")
    copyfile(bam, tmp_bam)
    copyfile(tsv, tmp_tsv)
    result = mhapi.type(tmp_bam, tmp_tsv, minbasequal=13)
    ac30 = result.data["markers"]["MHDBL000030"]["typing_result"]
    ac197 = result.data["markers"]["MHDBL000197"]["typing_result"]
    assert ac30 == {"A,A,T,C": 3, "A,C,C,C": 2, "A,C,C,G": 18, "G,C,C,C": 1, "G,C,C,G": 34}
    assert ac197 == {"A,A,T,T,T": 30, "A,A,T,T,C": 39, "A,A,T,A,T": 1, "A,A,T,A,C": 1}


def test_type_cli_simple(tmp_path):
    outfile = str(tmp_path / "typing-result.json")
    arglist = [
        "type",
        "--out",
        outfile,
        data_file("pashtun-sim/tiny-panel.tsv"),
        data_file("pashtun-sim/aligned-reads.bam"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.type.main(args)
    result = TypingResult(fromfile=outfile)
    assert result.haplotypes("mh13KK-218") == set()
    assert result.data["markers"]["mh13KK-218"]["typing_result"] == {
        "C,T,C,G": 1,
        "C,T,T,T": 1,
        "G,T,C,T": 1,
        "T,A,C,T": 1,
        "T,G,C,T": 3,
        "T,G,T,T": 2,
        "T,T,A,T": 5,
        "T,T,C,A": 1,
        "T,T,C,C": 2,
        "T,T,C,G": 2,
        "T,T,C,T": 1178,
        "T,T,G,T": 2,
        "T,T,T,A": 2,
        "T,T,T,G": 6,
        "T,T,T,T": 1170,
    }
    assert result.haplotypes("mh21KK-320") == set()
    assert result.data["markers"]["mh21KK-320"]["typing_result"] == {
        "G,A,A,A": 1,
        "G,A,C,A": 3,
        "G,A,G,A": 3,
        "G,A,T,A": 1075,
        "G,A,T,C": 1,
        "G,A,T,G": 1,
        "G,A,T,T": 2,
        "G,C,C,A": 1,
        "G,C,T,A": 4,
        "G,G,A,A": 2,
        "G,G,A,T": 1,
        "G,G,C,A": 1075,
        "G,G,C,C": 3,
        "G,G,C,G": 12,
        "G,G,C,T": 5,
        "G,G,T,A": 4,
        "G,T,C,A": 1,
        "T,G,C,A": 1,
    }


def test_type_filter_threshold():
    bam = data_file("bam/dyncut-test-reads.bam")
    tsv = data_file("def/dyncut-panel.tsv")
    rslt = mhapi.type(bam, tsv)
    rslt.filter(static=10, dynamic=0.005)
    assert rslt.haplotypes("MHDBL000018") == set(["C,A,C,T,G", "T,G,C,T,G"])
    assert rslt.haplotypes("MHDBL000156") == set(["T,C,A,C", "T,C,G,G"])
    rslt = mhapi.type(bam, tsv)
    rslt.filter(static=4, dynamic=0.005)
    assert rslt.haplotypes("MHDBL000018") == set(
        ["C,A,C,T,G", "T,G,C,T,G", "C,A,C,T,A", "T,G,C,T,A"]
    )
    assert rslt.haplotypes("MHDBL000156") == set(["T,C,A,C", "T,C,G,G"])


def test_type_no_var_offsets():
    bam = data_file("bam/sandawe-dad.bam")
    tsv = data_file("def/sandawe-empty.tsv")
    message = r"marker IDs unique to set1={mh01KK-205, mh02KK-005, mh03KK-006};"
    with pytest.raises(ValueError, match=message):
        mhapi.type(bam, tsv)


def test_type_marker_with_no_coverage():
    bam = data_file("bam/nocov.bam")
    tsv = data_file("def/nocov-offsets.tsv")
    result = mhapi.type(bam, tsv)
    for marker, mdata in result.data["markers"].items():
        assert "genotype" in mdata
