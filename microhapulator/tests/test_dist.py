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
from microhapulator.profile import TypingResult, SimulatedProfile
from microhapulator.tests import data_file
import pytest
from tempfile import NamedTemporaryFile


@pytest.mark.parametrize(
    "gt1,gt2,dist",
    [
        ("prof/gujarati-ind1-gt.json", "prof/gujarati-ind1-gt.json", 0),
        ("prof/gujarati-ind1-gt.json", "prof/gujarati-ind2-gt.json", 2),
        ("prof/gujarati-ind1-gt.json", "prof/gujarati-ind3-gt.json", 1),
        ("prof/gujarati-ind1-gt.json", "prof/gujarati-ind4-gt.json", 3),
        ("prof/gujarati-ind2-gt.json", "prof/gujarati-ind3-gt.json", 3),
        ("prof/gujarati-ind2-gt.json", "prof/gujarati-ind4-gt.json", 3),
        ("prof/gujarati-ind3-gt.json", "prof/gujarati-ind4-gt.json", 2),
    ],
)
def test_dist_gujarati(gt1, gt2, dist):
    r1 = TypingResult(data_file(gt1))
    r2 = TypingResult(data_file(gt2))
    assert microhapulator.op.dist(r1, r2) == dist


def test_dist_log_mixture():
    p1 = TypingResult(data_file("murica/y-obs-genotype.json"))
    p2 = SimulatedProfile.populate_from_bed(data_file("murica/y-sim-genotype.bed"))
    assert microhapulator.op.dist(p1, p2) == 19
    assert p1 != p2


def test_dist_even_mixture():
    with microhapulator.open(data_file("murica/x-obs-genotype.json"), "r") as fh:
        p1 = TypingResult(fh)
    p2 = SimulatedProfile.populate_from_bed(data_file("murica/x-sim-genotype.bed"))
    assert microhapulator.op.dist(p1, p2) == 0
    assert p1 == p2


@pytest.mark.parametrize("hdist", [0, 1, 2])
def test_dist_sim_vs_obs(hdist):
    with NamedTemporaryFile() as outfile:
        filename = "murica/z-obs-genotype-dist{:d}.json".format(hdist)
        arglist = [
            "dist",
            "--out",
            outfile.name,
            data_file(filename),
            data_file("murica/z-sim-genotype.json"),
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.cli.dist.main(args)
        with open(outfile.name, "r") as fh:
            assert json.load(fh) == {"hamming_distance": hdist}


def test_dist_cli():
    with NamedTemporaryFile() as outfile:
        arglist = [
            "dist",
            "--out",
            outfile.name,
            data_file("prof/gujarati-ind2-gt.json"),
            data_file("prof/gujarati-ind3-gt.json"),
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.cli.dist.main(args)
        with open(outfile.name, "r") as fh:
            assert json.load(fh) == {"hamming_distance": 3}
