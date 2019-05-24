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


@pytest.mark.parametrize('gt1,gt2,dist', [
    ('gujarati-ind1-gt.json', 'gujarati-ind1-gt.json', 0),
    ('gujarati-ind1-gt.json', 'gujarati-ind2-gt.json', 2),
    ('gujarati-ind1-gt.json', 'gujarati-ind3-gt.json', 1),
    ('gujarati-ind1-gt.json', 'gujarati-ind4-gt.json', 3),
    ('gujarati-ind2-gt.json', 'gujarati-ind3-gt.json', 3),
    ('gujarati-ind2-gt.json', 'gujarati-ind4-gt.json', 3),
    ('gujarati-ind3-gt.json', 'gujarati-ind4-gt.json', 2),
])
def test_dist_gujarati(gt1, gt2, dist):
    g1 = microhapulator.genotype.ObservedGenotype(data_file(gt1))
    g2 = microhapulator.genotype.ObservedGenotype(data_file(gt2))
    assert microhapulator.dist.dist(g1, g2) == dist


def test_dist_cli():
    with NamedTemporaryFile() as outfile:
        arglist = [
            'dist', '--out', outfile.name, '--json', data_file('gujarati-ind2-gt.json'),
            data_file('gujarati-ind3-gt.json')
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.dist.main(args)
        with open(outfile.name, 'r') as fh:
            assert json.load(fh) == {"hamming_distance": 3}
