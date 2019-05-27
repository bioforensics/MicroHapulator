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


def test_dist_log_mixture():
    f1 = data_file('murica/y-obs-genotype.json')
    g1 = microhapulator.genotype.ObservedGenotype(f1)
    f2 = data_file('murica/y-sim-genotype.bed')
    with microhapulator.open(f2, 'r') as fh:
        g2 = microhapulator.genotype.SimulatedGenotype(frombed=fh, ploidy=6)
    assert microhapulator.dist.dist(g1, g2) == 19
    assert g1 != g2


def test_dist_even_mixture():
    f1 = data_file('murica/x-obs-genotype.json')
    g1 = microhapulator.genotype.ObservedGenotype(f1)
    f2 = data_file('murica/x-sim-genotype.bed')
    with microhapulator.open(f2, 'r') as fh:
        g2 = microhapulator.genotype.SimulatedGenotype(frombed=fh, ploidy=6)
    assert microhapulator.dist.dist(g1, g2) == 0
    assert g1 == g2


@pytest.mark.parametrize('hdist', [0, 1, 2])
def test_dist_bed_vs_json(hdist):
    with NamedTemporaryFile() as outfile:
        filename = 'murica/z-obs-genotype-dist{:d}.json'.format(hdist)
        arglist = [
            'dist', '--out', outfile.name, '--obs', data_file(filename),
            '--sim', data_file('murica/z-sim-genotype.bed')
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.dist.main(args)
        with open(outfile.name, 'r') as fh:
            assert json.load(fh) == {"hamming_distance": hdist}


def test_dist_cli():
    with NamedTemporaryFile() as outfile:
        arglist = [
            'dist', '--out', outfile.name, '--obs', data_file('gujarati-ind2-gt.json'),
            data_file('gujarati-ind3-gt.json')
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.dist.main(args)
        with open(outfile.name, 'r') as fh:
            assert json.load(fh) == {"hamming_distance": 3}


@pytest.mark.parametrize('nargs', [1, 3, 4])
def test_dist_cli_bad(nargs):
    arglist = ['dist', '--obs'] + [data_file('gujarati-ind2-gt.json')] * nargs
    args = microhapulator.cli.parse_args(arglist)
    with pytest.raises(ValueError) as ve:
        microhapulator.dist.main(args)
    assert 'the number of BED and/or JSON files must be 2!' in str(ve)
