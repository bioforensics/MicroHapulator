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

import json
import microhapulator
from microhapulator.profile import ObservedProfile
from microhapulator.tests import data_file
import pytest
from shutil import copyfile
from tempfile import NamedTemporaryFile


def test_type_simple():
    bam = data_file("pashtun-sim/aligned-reads.bam")
    tsv = data_file("pashtun-sim/tiny-panel.tsv")
    observed = microhapulator.type.type(bam, tsv, static=10, dynamic=0.25)
    expected = ObservedProfile(fromfile=data_file("pashtun-sim/test-output.json"))
    assert observed == expected


def test_type_simpler():
    bam = data_file("pashtun-sim/aligned-reads.bam")
    tsv = data_file("pashtun-sim/tiny-panel.tsv")
    observed = microhapulator.type.type(bam, tsv)
    expected = ObservedProfile(fromfile=data_file("pashtun-sim/test-output-sans-genotype.json"))
    assert observed == expected


def test_type_missing_bam_index(tmp_path):
    bam = data_file("bam/three-contrib-log-link.bam")
    tsv = data_file("def/default-panel-offsets.tsv")
    tmp_bam = str(tmp_path / "reads.bam")
    tmp_tsv = str(tmp_path / "offsets.tsv")
    copyfile(bam, tmp_bam)
    copyfile(tsv, tmp_tsv)
    result = microhapulator.type.type(tmp_bam, tmp_tsv, minbasequal=13)
    ac30 = result.data["markers"]["MHDBL000030"]["allele_counts"]
    ac197 = result.data["markers"]["MHDBL000197"]["allele_counts"]
    assert ac30 == {"A,A,T,C": 3, "A,C,C,C": 2, "A,C,C,G": 18, "G,C,C,C": 1, "G,C,C,G": 34}
    assert ac197 == {"A,A,T,T,T": 30, "A,A,T,T,C": 39, "A,A,T,A,T": 1, "A,A,T,A,C": 1}


def test_type_cli_simple(tmp_path):
    outfile = str(tmp_path / "typing-result.json")
    arglist = [
        "type",
        "--out",
        outfile,
        "--static",
        "5",
        "--dynamic",
        "0.25",
        data_file("pashtun-sim/tiny-panel.tsv"),
        data_file("pashtun-sim/aligned-reads.bam"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.type.main(args)
    observed = ObservedProfile(fromfile=outfile)
    expected = ObservedProfile(fromfile=data_file("pashtun-sim/test-output.json"))
    assert observed == expected


def test_type_dyn_cutoff():
    bam = data_file("bam/dyncut-test-reads.bam")
    tsv = data_file("def/dyncut-panel.tsv")
    rslt = microhapulator.type.type(bam, tsv, static=10, dynamic=0.25)
    assert rslt.alleles("MHDBL000018") == set(["C,A,C,T,G", "T,G,C,T,G"])
    assert rslt.alleles("MHDBL000156") == set(["T,C,A,C", "T,C,G,G"])
    rslt = microhapulator.type.type(bam, tsv, static=4, dynamic=0.25)
    assert rslt.alleles("MHDBL000018") == set(["C,A,C,T,G", "T,G,C,T,G", "C,A,C,T,A", "T,G,C,T,A"])
    assert rslt.alleles("MHDBL000156") == set(["T,C,A,C", "T,C,G,G"])


def test_type_no_var_offsets():
    bam = data_file("bam/sandawe-dad.bam")
    tsv = data_file("def/sandawe-empty.tsv")
    message = r"marker IDs unique to set1={mh01KK-205, mh02KK-005, mh03KK-006};"
    with pytest.raises(ValueError, match=message):
        result = microhapulator.type.type(bam, tsv)
