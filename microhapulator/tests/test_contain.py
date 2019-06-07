#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
import microhapulator
from microhapulator.tests import data_file
import pytest
from tempfile import NamedTemporaryFile


@pytest.mark.parametrize('f1,f2,total,contained', [
    ('four-brits-sim.json', 'one-brit-sim.json', 38, 38),
    ('four-brits-sim.json', 'one-american-sim.json', 38, 31),
    ('one-brit-sim.json', 'one-american-sim.json', 38, 11),
])
def test_contain(f1, f2, total, contained):
    gt1 = microhapulator.genotype.Genotype(data_file(f1))
    gt2 = microhapulator.genotype.Genotype(data_file(f2))
    c, t = microhapulator.contain.contain(gt1, gt2)
    assert t == total
    assert c == contained


def test_contain_cli(capsys):
    arglist = [
        'contain', data_file('one-brit-sim.json'), data_file('one-italian-sim.json')
    ]
    args = microhapulator.cli.parse_args(arglist)
    microhapulator.contain.main(args)
    terminal = capsys.readouterr()
    assert '"containment": 0.4444' in terminal.out
