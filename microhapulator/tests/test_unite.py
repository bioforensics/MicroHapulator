#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import numpy.random
import microhapulator
from microhapulator.genotype import Genotype
from microhapulator.tests import data_file
import pytest


@pytest.mark.parametrize('momgt,dadgt,kidgt,seed', [
    ('green-mom-1-gt.json', 'green-dad-1-gt.json', 'green-kid-1-gt.json', 110517),
    ('green-mom-2-gt.json', 'green-dad-2-gt.json', 'green-kid-2-gt.json', 111017),
])
def test_unite_basic(momgt, dadgt, kidgt, seed):
    mom = Genotype(fromfile=data_file(momgt))
    dad = Genotype(fromfile=data_file(dadgt))
    kid = Genotype(fromfile=data_file(kidgt))
    numpy.random.seed(seed)
    test = Genotype.unite(mom, dad)
    test.dump('BOOGER' + kidgt)
    assert test == kid
