#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
from microhapulator.tests import data_file
import pytest


def test_rmp():
    gt = microhapulator.genotype.Genotype(fromfile=data_file('korea-5loc.json'))
    assert gt.rand_match_prob('MHDBP000053') == pytest.approx(9.30529727667e-10)


def test_rmp_lrt():
    gt1 = microhapulator.genotype.Genotype(fromfile=data_file('korea-5loc.json'))
    gt2 = microhapulator.genotype.Genotype(fromfile=data_file('korea-5loc-1diff.json'))
    assert gt1.rmp_lr_test(gt1, 'MHDBP000053') == pytest.approx(1074656693.1355)
    assert gt1.rmp_lr_test(gt2, 'MHDBP000053') == pytest.approx(1074656.6931355)


@pytest.mark.parametrize('altfile,lrvalue', [
    ('korea-5loc-2diff-a.json', 1943.0058),
    ('korea-5loc-2diff-b.json', 9624.2067),
    ('korea-5loc-2diff-c.json', 34191.2606),
])
def test_rmp_lrt_2diff(altfile, lrvalue):
    gt1 = microhapulator.genotype.Genotype(fromfile=data_file('korea-5loc.json'))
    gt2 = microhapulator.genotype.Genotype(fromfile=data_file(altfile))
    assert gt1.rmp_lr_test(gt2, 'MHDBP000053') == pytest.approx(1074.6567)
    assert gt2.rmp_lr_test(gt1, 'MHDBP000053') == pytest.approx(lrvalue)


def test_prob_cli_rmp(capsys):
    arglist = [
        'prob', 'MHDBP000053', data_file('korea-5loc.json')
    ]
    args = microhapulator.cli.parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    assert '"random_match_probability": "9.305E-10"' in terminal.out


def test_prob_cli_lrt(capsys):
    arglist = [
        'prob', 'MHDBP000053', data_file('korea-5loc.json'), data_file('korea-5loc-1diff.json')
    ]
    args = microhapulator.cli.parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    assert '"rmp_likelihood_ratio": "1.075E+06"' in terminal.out


def test_prob_zero_freq():
    gt = microhapulator.genotype.Genotype(fromfile=data_file('korea-5loc-zerofreq.json'))
    assert gt.rand_match_prob('MHDBP000053') == pytest.approx(2.963E-12)


def test_prob_missing_freq():
    gt = microhapulator.genotype.Genotype(fromfile=data_file('korea-5loc-missfreq.json'))
    assert gt.rand_match_prob('MHDBP000053') == pytest.approx(4.898E-11)
