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
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import pandas as pd
import pytest


@pytest.fixture(scope="session")
def k5freqs():
    return pd.read_csv(data_file("freq/korea-5loc-freq.tsv"), sep="\t")


def test_rmp(k5freqs):
    p = Profile(fromfile=data_file("prof/korea-5loc.json"))
    assert p.rand_match_prob(k5freqs) == pytest.approx(7.444e-09)


def test_rmp_lrt(k5freqs):
    p1 = Profile(fromfile=data_file("prof/korea-5loc.json"))
    p2 = Profile(fromfile=data_file("prof/korea-5loc-1diff.json"))
    assert p1.rmp_lr_test(p1, k5freqs) == pytest.approx(134332086.64194357)
    assert p1.rmp_lr_test(p2, k5freqs) == pytest.approx(134332.08664194357)


@pytest.mark.parametrize(
    "altfile,lrvalue",
    [
        ("prof/korea-5loc-2diff-a.json", 121.4379),
        ("prof/korea-5loc-2diff-b.json", 1203.0258),
        ("prof/korea-5loc-2diff-c.json", 2136.9538),
    ],
)
def test_rmp_lrt_2diff(altfile, lrvalue, k5freqs):
    p1 = Profile(fromfile=data_file("prof/korea-5loc.json"))
    p2 = Profile(fromfile=data_file(altfile))
    assert p1.rmp_lr_test(p2, k5freqs) == pytest.approx(134.3321)
    assert p2.rmp_lr_test(p1, k5freqs) == pytest.approx(lrvalue)


def test_prob_cli_rmp(capsys):
    arglist = ["prob", data_file("freq/korea-5loc-freq.tsv"), data_file("prof/korea-5loc.json")]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert '"random_match_probability": "7.444E-09"' in terminal.out


def test_prob_cli_lrt(capsys):
    arglist = [
        "prob",
        data_file("freq/korea-5loc-freq.tsv"),
        data_file("prof/korea-5loc.json"),
        data_file("prof/korea-5loc-1diff.json"),
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert '"likelihood_ratio": "1.343E+05"' in terminal.out


def test_prob_zero_freq(k5freqs):
    p = Profile(fromfile=data_file("prof/korea-5loc-zerofreq.json"))
    assert p.rand_match_prob(k5freqs) == pytest.approx(2.3708e-11)


def test_prob_missing_freq(k5freqs):
    p = Profile(fromfile=data_file("prof/korea-5loc-missfreq.json"))
    assert p.rand_match_prob(k5freqs) == pytest.approx(7.8360e-10)
