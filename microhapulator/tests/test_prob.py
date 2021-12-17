#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapdb
import microhapulator
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import pytest

FREQS = microhapdb.frequencies[microhapdb.frequencies.Population == "SA000936S"]


def test_rmp():
    p = Profile(fromfile=data_file("korea-5loc.json"))
    assert p.rand_match_prob(FREQS) == pytest.approx(7.444e-09)


def test_rmp_lrt():
    p1 = Profile(fromfile=data_file("korea-5loc.json"))
    p2 = Profile(fromfile=data_file("korea-5loc-1diff.json"))
    assert p1.rmp_lr_test(p1, FREQS) == pytest.approx(134332086.64194357)
    assert p1.rmp_lr_test(p2, FREQS) == pytest.approx(134332.08664194357)


@pytest.mark.parametrize(
    "altfile,lrvalue",
    [
        ("korea-5loc-2diff-a.json", 121.4379),
        ("korea-5loc-2diff-b.json", 1203.0258),
        ("korea-5loc-2diff-c.json", 2136.9538),
    ],
)
def test_rmp_lrt_2diff(altfile, lrvalue):
    p1 = Profile(fromfile=data_file("korea-5loc.json"))
    p2 = Profile(fromfile=data_file(altfile))
    assert p1.rmp_lr_test(p2, FREQS) == pytest.approx(134.3321)
    assert p2.rmp_lr_test(p1, FREQS) == pytest.approx(lrvalue)


def test_prob_cli_rmp(capsys):
    arglist = ["prob", "SA000936S", data_file("korea-5loc.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert '"random_match_probability": "7.444E-09"' in terminal.out


def test_prob_cli_lrt(capsys):
    arglist = [
        "prob",
        "SA000936S",
        data_file("korea-5loc.json"),
        data_file("korea-5loc-1diff.json"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert '"likelihood_ratio": "1.343E+05"' in terminal.out


def test_prob_zero_freq():
    p = Profile(fromfile=data_file("korea-5loc-zerofreq.json"))
    assert p.rand_match_prob(FREQS) == pytest.approx(2.3708e-11)


def test_prob_missing_freq():
    p = Profile(fromfile=data_file("korea-5loc-missfreq.json"))
    assert p.rand_match_prob(FREQS) == pytest.approx(7.8360e-10)


def test_bad_pop():
    arglist = ["prob", "Han", data_file("korea-5loc.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    message = r'issue with population "Han"; invalid or not unique'
    with pytest.raises(ValueError, match=message):
        microhapulator.prob.main(args)

    arglist = ["prob", "FakePopulation", data_file("korea-5loc.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    message = r'issue with population "FakePopulation"; invalid or not unique'
    with pytest.raises(ValueError, match=message):
        microhapulator.prob.main(args)
