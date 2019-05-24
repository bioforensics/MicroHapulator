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
from microhapulator.type import MissingBAMIndexError
import pytest
from tempfile import NamedTemporaryFile


def test_type_simple():
    bam = data_file('pashtun-sim/aligned-reads.bam')
    fasta = data_file('pashtun-sim/tiny-panel.fasta.gz')
    gt = microhapulator.type.type(bam, fasta)
    testgtfile = data_file('pashtun-sim/test-output.json')
    testgt = microhapulator.genotype.ObservedGenotype(filename=testgtfile)
    assert gt.data == testgt.data
    assert gt.dump() == testgt.dump()


def test_type_missing_bam_index():
    bam = data_file('three-contrib-log-link.bam')
    fasta = data_file('default-panel.fasta.gz')
    with pytest.raises(MissingBAMIndexError) as ie:
        gt = microhapulator.type.type(bam, fasta)
    assert 'Please index' in str(ie)


def test_type_cli_simple():
    with NamedTemporaryFile() as outfile:
        arglist = [
            'type', '--out', outfile.name, data_file('pashtun-sim/tiny-panel.fasta.gz'),
            data_file('pashtun-sim/aligned-reads.bam'),
        ]
        args = microhapulator.cli.parse_args(arglist)
        microhapulator.type.main(args)
        testgtfile = data_file('pashtun-sim/test-output.json')
        with open(outfile.name, 'r') as fh, open(testgtfile, 'r') as testfh:
            gtdata = json.load(fh)
            testgtdata = json.load(testfh)
            assert gtdata == testgtdata
