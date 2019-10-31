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
from microhapulator.profile import Profile
from microhapulator.tests import data_file
import pytest
from tempfile import NamedTemporaryFile


@pytest.mark.parametrize('momgt,dadgt,kidgt,seed', [
    ('green-mom-1-gt.json', 'green-dad-1-gt.json', 'green-kid-1-gt.json', 110517),
    ('green-mom-2-gt.json', 'green-dad-2-gt.json', 'green-kid-2-gt.json', 111017),
])
def test_unite_basic(momgt, dadgt, kidgt, seed):
    mom = Profile(fromfile=data_file(momgt))
    dad = Profile(fromfile=data_file(dadgt))
    kid = Profile(fromfile=data_file(kidgt))
    numpy.random.seed(seed)
    test = Profile.unite(mom, dad)
    assert test == kid


def test_unite_cli():
    with NamedTemporaryFile(suffix='.json') as outfile:
        arglist = [
            'unite', '--seed', '113817', '--out', outfile.name,
            data_file('green-mom-3-gt.json'), data_file('green-dad-3-gt.json'),
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.unite.main(args)
        p = Profile(fromfile=outfile.name)
        testp = Profile(fromfile=data_file('green-kid-3-gt.json'))
        assert p == testp
