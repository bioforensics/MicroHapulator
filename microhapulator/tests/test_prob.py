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
    assert gt.rand_match_prob('MHDBP000053') == pytest.approx(1.2122586342e-07)


def test_rmp_lrt():
    gt1 = microhapulator.genotype.Genotype(fromfile=data_file('korea-5loc.json'))
    gt2 = microhapulator.genotype.Genotype(fromfile=data_file('korea-5loc-1diff.json'))
    assert gt1.rmp_lr_test(gt1, 'MHDBP000053') == pytest.approx(8249064.776508)
    assert gt1.rmp_lr_test(gt2, 'MHDBP000053') == pytest.approx(7424158.298857)


def test_prob_cli_rmp(capsys):
    arglist = [
        'prob', 'MHDBP000053', data_file('korea-5loc.json')
    ]
    args = microhapulator.cli.parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    assert '"random_match_probability": "1.212E-07"' in terminal.out


def test_prob_cli_lrt(capsys):
    arglist = [
        'prob', 'MHDBP000053', data_file('korea-5loc.json'), data_file('korea-5loc-1diff.json')
    ]
    args = microhapulator.cli.parse_args(arglist)
    microhapulator.prob.main(args)
    terminal = capsys.readouterr()
    assert '"rmp_likelihood_ratio": "7.424E+06"' in terminal.out
