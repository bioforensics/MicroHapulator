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
from microhapulator.thresholds import ThresholdIndex
import pytest
from shutil import copyfile


def test_type_simple():
    bam = data_file("pashtun-sim/aligned-reads.bam")
    tsv = data_file("pashtun-sim/tiny-panel.tsv")
    result = mhapi.type(bam, tsv)
    assert result.haplotypes("mh13KK-218") == set()
    assert result.data["markers"]["mh13KK-218"]["typing_result"] == {
        "T,T,C,T": 1178,
        "T,T,T,T": 1170,
        "T,T,T,G": 6,
        "T,T,C,G": 2,
        "G,T,C,T": 1,
        "T,T,C,A": 1,
        "T,G,T,T": 2,
        "T,T,A,T": 5,
        "T,G,C,T": 3,
        "T,A,C,T": 1,
        "C,T,C,G": 1,
        "T,T,G,T": 2,
        "T,T,C,C": 2,
        "T,T,T,A": 2,
        "T,T,C,-": 1,
        "C,T,T,T": 1,
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
        "G,A,G,-": 1,
        "G,G,A,T": 1,
        "G,G,C,A": 1075,
        "G,G,C,C": 3,
        "G,G,C,G": 12,
        "G,G,C,T": 5,
        "G,G,T,A": 4,
        "G,T,C,A": 1,
        "T,G,C,A": 1,
    }


def test_type_multiple_marker_definitions_per_locus():
    bam = data_file("pashtun-sim/aligned-reads-multi.bam")
    tsv = data_file("pashtun-sim/tiny-panel-multidef.tsv")
    result = mhapi.type(bam, tsv)
    assert sorted(result.data["markers"]) == [
        "mh13KK-218.v1",
        "mh13KK-218.v6",
        "mh21KK-320.v1",
        "mh21KK-320.v5",
    ]
    assert result.data["markers"]["mh21KK-320.v1"]["typing_result"]["G,A,T,A"] == 1075
    assert result.data["markers"]["mh21KK-320.v1"]["typing_result"]["G,A,C,A"] == 3
    assert result.data["markers"]["mh21KK-320.v5"]["typing_result"]["G,A,T,A,A"] == 932
    assert result.data["markers"]["mh21KK-320.v5"]["typing_result"]["G,G,A,T,A"] == 1


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
        "T,T,C,-": 1,
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
        "G,A,G,-": 1,
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


@pytest.mark.parametrize(
    "static,dynamic,genotype",
    [
        (5, 0.00005, {"C,A,G,A", "C,G,G,A", "C,T,G,A"}),
        (10, 0.00005, {"C,A,G,A", "C,G,G,A"}),
        (5, 0.02, {"C,A,G,A", "C,G,G,A"}),
        (2000, 0.02, {"C,A,G,A"}),
    ],
)
def test_type_filter_threshold(static, dynamic, genotype):
    result = TypingResult(fromfile=data_file("prof/srm.json"))
    thresholds = ThresholdIndex.load(global_static=static, global_dynamic=dynamic)
    result.filter(thresholds)
    assert result.haplotypes("mh05KK-170.v1") == genotype


def test_type_no_var_offsets():
    bam = data_file("bam/sandawe-dad.bam")
    tsv = data_file("def/sandawe-empty.tsv")
    message = "no marker definitions in microhap index"
    with pytest.warns(UserWarning, match=message):
        mhapi.type(bam, tsv)


def test_type_marker_with_no_coverage():
    bam = data_file("bam/nocov.bam")
    tsv = data_file("def/nocov-offsets.tsv")
    result = mhapi.type(bam, tsv)
    for marker, mdata in result.data["markers"].items():
        assert "genotype" in mdata


def test_type_marker_has_indel_spanning_snp():
    bam = data_file("bam/indel_snp.bam")
    tsv = data_file("def/indel_snp.tsv")
    result = mhapi.type(bam, tsv)
    assert "mh09WL-026" in result.data["markers"]
    for allele in result.data["markers"]["mh09WL-026"]["typing_result"]:
        snp_alleles = allele.split(",")
        assert len(snp_alleles) == 6
        assert snp_alleles[-2] == "-"
