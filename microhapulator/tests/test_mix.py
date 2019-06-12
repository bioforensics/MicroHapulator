#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
from microhapulator.genotype import SimulatedGenotype
from microhapulator.tests import data_file
import pytest
from tempfile import NamedTemporaryFile


def test_mix_main():
    with NamedTemporaryFile(suffix='.json.gz') as outfile:
        arglist = [
            'mix', '--out', outfile.name, data_file('green-sim-gt-1.json.gz'),
            data_file('green-sim-gt-2.json.gz'), data_file('green-sim-gt-3.json.gz')
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.mix.main(args)
        gt = SimulatedGenotype(fromfile=outfile.name)
        testgt = SimulatedGenotype(fromfile=data_file('green-sim-gt-combined.json.gz'))
        assert gt == testgt
