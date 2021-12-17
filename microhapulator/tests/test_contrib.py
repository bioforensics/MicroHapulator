#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import pytest


@pytest.mark.parametrize(
    "pjson,numcontrib",
    [
        ("single-contrib-1.json", 1),
        ("single-contrib-2.json", 1),
        ("single-contrib-3.json", 1),
        ("two-contrib-even.json", 2),
        ("three-contrib-even.json", 3),
        ("three-contrib-log.json", 3),
    ],
)
def test_contrib_json(pjson, numcontrib):
    profile = microhapulator.contrib.load_profile(json=data_file(pjson))
    n, *data = microhapulator.contrib.contrib(profile)
    assert n == numcontrib


def test_contrib_bam():
    bam = data_file("three-contrib-log.bam")
    refr = data_file("default-panel.fasta.gz")
    profile = microhapulator.contrib.load_profile(
        bamfile=bam, refrfasta=refr, dynamic=0.25, static=10
    )
    n, *data = microhapulator.contrib.contrib(profile)
    assert n == 3


def test_contrib_main(capsys):
    bam = data_file("three-contrib-log.bam")
    refr = data_file("default-panel.fasta.gz")
    arglist = ["contrib", "-b", bam, "-r", refr, "--static", "10", "--dynamic", "0.25"]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.contrib.main(args)
    out, err = capsys.readouterr()
    assert '"min_num_contrib": 3' in out


def test_main_no_op():
    arglist = ["contrib"]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    pattern = r"must provide either JSON profile or BAM and refr FASTA"
    with pytest.raises(ValueError, match=pattern) as ve:
        microhapulator.contrib.main(args)
