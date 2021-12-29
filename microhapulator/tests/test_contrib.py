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
from microhapulator.tests import data_file
import pytest


@pytest.mark.parametrize(
    "pjson,numcontrib",
    [
        ("prof/single-contrib-1.json", 1),
        ("prof/single-contrib-2.json", 1),
        ("prof/single-contrib-3.json", 1),
        ("prof/two-contrib-even.json", 2),
        ("prof/three-contrib-even.json", 3),
        ("prof/three-contrib-log.json", 3),
    ],
)
def test_contrib_json(pjson, numcontrib):
    profile = microhapulator.cli.contrib.load_profile(json=data_file(pjson))
    n, *data = mhapi.contrib(profile)
    assert n == numcontrib


def test_contrib_bam():
    bam = data_file("bam/three-contrib-log.bam")
    defn = data_file("def/default-panel-offsets.tsv")
    profile = microhapulator.cli.contrib.load_profile(
        bamfile=bam, markertsv=defn, dynamic=0.25, static=10
    )
    n, *data = mhapi.contrib(profile)
    assert n == 3


def test_contrib_main(capsys):
    bam = data_file("bam/three-contrib-log.bam")
    defn = data_file("def/default-panel-offsets.tsv")
    arglist = ["contrib", "-b", bam, "-t", defn, "--static", "10", "--dynamic", "0.25"]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.cli.contrib.main(args)
    out, err = capsys.readouterr()
    assert '"min_num_contrib": 3' in out


def test_main_no_op():
    arglist = ["contrib"]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    pattern = r"must provide either JSON profile or BAM and refr FASTA"
    with pytest.raises(ValueError, match=pattern) as ve:
        microhapulator.cli.contrib.main(args)
