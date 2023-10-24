# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
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
from microhapulator import load_marker_thresholds
import microhapulator.api as mhapi
from microhapulator.profile import TypingResult
from microhapulator.tests import data_file
import pandas as pd
import pytest


def test_filter_simple():
    bam = data_file("pashtun-sim/aligned-reads.bam")
    tsv = data_file("pashtun-sim/tiny-panel.tsv")
    observed = mhapi.type(bam, tsv)
    thresholds = load_marker_thresholds(observed.markers(), global_static=10, global_dynamic=0.05)
    observed.filter(thresholds)
    expected = TypingResult(fromfile=data_file("pashtun-sim/test-output.json"))
    assert observed == expected


def test_filter_config_file():
    result = TypingResult(fromfile=data_file("prof/deep-raw.json"))
    thresholds = load_marker_thresholds(
        result.markers(), configfile=data_file("filters.csv"), global_static=5, global_dynamic=0.02
    )
    result.filter(thresholds)
    assert len(result.haplotypes("mh01XYZ-1")) == 8
    assert len(result.haplotypes("mh02XYZ-2")) == 2
    assert len(result.haplotypes("mh02XYZ-3")) == 2


def test_filter_missing_column():
    result = TypingResult(fromfile=data_file("prof/deep-raw.json"))
    with pytest.raises(ValueError, match=r"filter config file missing column\(s\): Static"):
        thresholds = load_marker_thresholds(
            result.markers(),
            configfile=data_file("filters-missing.csv"),
            global_static=5,
            global_dynamic=0.02,
        )


def test_filter_dupl_marker():
    result = TypingResult(fromfile=data_file("prof/deep-raw.json"))
    message = "filter config file contains duplicate entries for some markers"
    with pytest.raises(ValueError, match=message):
        thresholds = load_marker_thresholds(
            result.markers(),
            configfile=data_file("filters-redundant.csv"),
            global_static=5,
            global_dynamic=0.02,
        )


def test_filter_cli(tmp_path):
    unfiltered = str(tmp_path / "typing-result.json")
    filtered = str(tmp_path / "genotype-call.json")
    arglist = [
        "type",
        "--out",
        unfiltered,
        data_file("pashtun-sim/tiny-panel.tsv"),
        data_file("pashtun-sim/aligned-reads.bam"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.type.main(args)
    arglist = ["filter", "--out", filtered, "--static", "5", "--dynamic", "0.05", unfiltered]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.filter.main(args)
    observed = TypingResult(fromfile=filtered)
    expected = TypingResult(fromfile=data_file("pashtun-sim/test-output.json"))
    assert observed == expected


def test_filter_cli_config(tmp_path):
    unfiltered = data_file("prof/deep-raw.json")
    filtered = str(tmp_path / "genotype-call.json")
    arglist = [
        "filter",
        unfiltered,
        "--out",
        filtered,
        "--static",
        "5",
        "--dynamic",
        "0.02",
        "--config",
        data_file("filters.csv"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.filter.main(args)
    observed = TypingResult(fromfile=filtered)
    expected = TypingResult(fromfile=data_file("prof/deep-filt.json"))
    assert observed == expected
