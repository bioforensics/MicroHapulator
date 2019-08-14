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
    testgt = microhapulator.genotype.ObservedGenotype(fromfile=testgtfile)
    assert gt == testgt


def test_type_missing_bam_index():
    bam = data_file('three-contrib-log-link.bam')
    fasta = data_file('default-panel.fasta.gz')
    with pytest.raises(MissingBAMIndexError) as ie:
        gt = microhapulator.type.type(bam, fasta)
    assert 'Please index' in str(ie)


def test_type_cli_simple():
    with NamedTemporaryFile() as outfile:
        arglist = [
            'type', '--out', outfile.name, '--threshold', '5',
            data_file('pashtun-sim/tiny-panel.fasta.gz'),
            data_file('pashtun-sim/aligned-reads.bam'),
        ]
        args = microhapulator.cli.get_parser().parse_args(arglist)
        microhapulator.type.main(args)
        testgtfile = data_file('pashtun-sim/test-output.json')
        gtdata = microhapulator.genotype.ObservedGenotype(fromfile=outfile.name)
        testgtdata = microhapulator.genotype.ObservedGenotype(fromfile=testgtfile)
        assert gtdata == testgtdata


def test_type_dyn_cutoff():
    bam = data_file('dyncut-test-reads.bam')
    fasta = data_file('dyncut-panel.fasta.gz')

    gt = microhapulator.type.type(bam, fasta)
    assert gt.alleles('MHDBL000018') == set(['C,A,C,T,G', 'T,G,C,T,G'])
    assert gt.alleles('MHDBL000156') == set(['T,C,A,C', 'T,C,G,G'])

    gt = microhapulator.type.type(bam, fasta, threshold=4)
    assert gt.alleles('MHDBL000018') == set(['C,A,C,T,G', 'T,G,C,T,G', 'C,A,C,T,A', 'T,G,C,T,A'])
    assert gt.alleles('MHDBL000156') == set(['T,C,A,C', 'T,C,G,G'])
