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


def test_rmp():
    p = Profile(fromfile=data_file('korea-5loc.json'))
    assert p.rand_match_prob('SA000936S') == pytest.approx(9.30529727667e-10)


def test_rmp_lrt():
    p1 = Profile(fromfile=data_file('korea-5loc.json'))
    p2 = Profile(fromfile=data_file('korea-5loc-1diff.json'))
    assert p1.rmp_lr_test(p1, 'SA000936S') == pytest.approx(1074656693.1355)
    assert p1.rmp_lr_test(p2, 'SA000936S') == pytest.approx(1074656.6931355)


@pytest.mark.parametrize('altfile,lrvalue', [
    ('korea-5loc-2diff-a.json', 1943.0058),
    ('korea-5loc-2diff-b.json', 9624.2067),
    ('korea-5loc-2diff-c.json', 34191.2606),
])
def test_rmp_lrt_2diff(altfile, lrvalue):
    p1 = Profile(fromfile=data_file('korea-5loc.json'))
    p2 = Profile(fromfile=data_file(altfile))
    assert p1.rmp_lr_test(p2, 'SA000936S') == pytest.approx(1074.6567)
    assert p2.rmp_lr_test(p1, 'SA000936S') == pytest.approx(lrvalue)


def test_prob_cli_rmp(capsys):
    arglist = [
        'prob', 'SA000936S', data_file('korea-5loc.json')
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    assert '"random_match_probability": "9.305E-10"' in terminal.out


def test_prob_cli_lrt(capsys):
    arglist = [
        'prob', 'SA000936S', data_file('korea-5loc.json'), data_file('korea-5loc-1diff.json')
    ]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    assert '"likelihood_ratio": "1.075E+06"' in terminal.out


def test_prob_zero_freq():
    p = Profile(fromfile=data_file('korea-5loc-zerofreq.json'))
    assert p.rand_match_prob('SA000936S') == pytest.approx(2.963E-12)


def test_prob_missing_freq():
    p = Profile(fromfile=data_file('korea-5loc-missfreq.json'))
    assert p.rand_match_prob('SA000936S') == pytest.approx(4.898E-11)


def test_bad_pop():
    arglist = ['prob', 'Han', data_file('korea-5loc.json')]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    message = r'issue with population "Han"; invalid or not unique'
    with pytest.raises(ValueError, match=message):
        microhapulator.prob.main(args)

    arglist = ['prob', 'FakePopulation', data_file('korea-5loc.json')]
    args = microhapulator.cli.get_parser().parse_args(arglist)
    message = r'issue with population "FakePopulation"; invalid or not unique'
    with pytest.raises(ValueError, match=message):
        microhapulator.prob.main(args)
