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
from microhapulator.profile import ObservedProfile
from microhapulator.tests import data_file
from microhapulator.type import MissingBAMIndexError
import pytest
from tempfile import NamedTemporaryFile


def test_type_simple():
    bam = data_file('pashtun-sim/aligned-reads.bam')
    fasta = data_file('pashtun-sim/tiny-panel.fasta.gz')
    gt = microhapulator.type.type(bam, fasta)
    testgtfile = data_file('pashtun-sim/test-output.json')
    testgt = ObservedProfile(fromfile=testgtfile)
    assert gt == testgt


def test_type_missing_bam_index():
    bam = data_file('three-contrib-log-link.bam')
    fasta = data_file('default-panel.fasta.gz')
    with pytest.raises(MissingBAMIndexError, match=r'Please index') as ie:
        gt = microhapulator.type.type(bam, fasta)


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
        gtdata = ObservedProfile(fromfile=outfile.name)
        testgtdata = ObservedProfile(fromfile=testgtfile)
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


def test_type_no_var_offsets():
    bam = data_file('sandawe-dad.bam')
    fasta = data_file('sandawe-panel.fasta.gz')
    with pytest.raises(ValueError, match=r'variant offsets not annotated for target amplicon:'):
        profile = microhapulator.type.type(bam, fasta)
