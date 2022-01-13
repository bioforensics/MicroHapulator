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
import microhapulator.api as mhapi
from microhapulator.profile import TypingResult
from microhapulator.tests import data_file
import pytest


def test_filter_simple():
    bam = data_file("pashtun-sim/aligned-reads.bam")
    tsv = data_file("pashtun-sim/tiny-panel.tsv")
    observed = mhapi.type(bam, tsv)
    observed.infer(static=10, dynamic=0.25)
    expected = TypingResult(fromfile=data_file("pashtun-sim/test-output.json"))
    assert observed == expected


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
    arglist = ["filter", "--out", filtered, "--static", "5", "--dynamic", "0.25", unfiltered]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.filter.main(args)
    observed = TypingResult(fromfile=filtered)
    expected = TypingResult(fromfile=data_file("pashtun-sim/test-output.json"))
    assert observed == expected
