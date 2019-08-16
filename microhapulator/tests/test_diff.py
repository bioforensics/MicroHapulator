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


def test_diff_basic():
    gt1 = microhapulator.genotype.SimulatedGenotype(fromfile=data_file('diff-comp-1.json'))
    gt2 = microhapulator.genotype.SimulatedGenotype(fromfile=data_file('diff-comp-2.json'))
    diff = list(microhapulator.diff.diff(gt1, gt2))
    assert diff == [
        ('MHDBL000140', {'C,C,A,A'}, {'C,C,T,A'}),
        ('MHDBL000163', {'A,A,G,A,T'}, {'C,G,A,A,T'}),
    ]


def test_diff_large():
    gt1 = microhapulator.genotype.SimulatedGenotype(fromfile=data_file('diff-comp-1.json'))
    gt2 = microhapulator.genotype.SimulatedGenotype(fromfile=data_file('diff-comp-3.json'))
    diff = list(microhapulator.diff.diff(gt1, gt2))
    loci = [d[0] for d in diff]
    print(diff[9], diff[13], diff[16])
    assert loci == [
        'MHDBL000002', 'MHDBL000007', 'MHDBL000013', 'MHDBL000017', 'MHDBL000018', 'MHDBL000030',
        'MHDBL000036', 'MHDBL000038', 'MHDBL000058', 'MHDBL000076', 'MHDBL000079', 'MHDBL000082',
        'MHDBL000085', 'MHDBL000101', 'MHDBL000106', 'MHDBL000111', 'MHDBL000112', 'MHDBL000122',
        'MHDBL000128', 'MHDBL000129', 'MHDBL000135', 'MHDBL000138', 'MHDBL000140', 'MHDBL000144',
        'MHDBL000154', 'MHDBL000163', 'MHDBL000183', 'MHDBL000210', 'MHDBL000211', 'MHDBL000212'
    ]
    assert diff[9] == ('MHDBL000076', {'G,T'}, set())
    assert diff[13] == ('MHDBL000101', {'C,C,C,T'}, {'T,C,C,C'})
    assert diff[16] == ('MHDBL000112', {'G,G,A,C'}, set())


def test_diff_cli():
    f1 = data_file('diff-comp-1.json')
    f2 = data_file('diff-comp-3.json')
    with NamedTemporaryFile(suffix='.json') as outfile:
        arglist = ['diff', '-o', outfile.name, f1, f2]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.diff.main(args)
        with microhapulator.open(outfile.name, 'r') as fh:
            output = fh.read().strip()
        with microhapulator.open(data_file('diff-comp-1-3.txt'), 'r') as fh:
            testoutput = fh.read().strip()
        assert output == testoutput
