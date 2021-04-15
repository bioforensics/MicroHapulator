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
import pytest
from shutil import copyfile
from tempfile import NamedTemporaryFile


def test_type_simple():
    bam = data_file('pashtun-sim/aligned-reads.bam')
    fasta = data_file('pashtun-sim/tiny-panel.fasta.gz')
    gt = microhapulator.type.type(bam, fasta, static=10)
    testgtfile = data_file('pashtun-sim/test-output.json')
    testgt = ObservedProfile(fromfile=testgtfile)
    assert gt == testgt


def test_type_simpler():
    bam = data_file('pashtun-sim/aligned-reads.bam')
    fasta = data_file('pashtun-sim/tiny-panel.fasta.gz')
    gt = microhapulator.type.type(bam, fasta)
    testgtfile = data_file('pashtun-sim/test-output-sans-genotype.json')
    testgt = ObservedProfile(fromfile=testgtfile)
    assert gt == testgt


def test_type_missing_bam_index(tmp_path):
    bam = data_file('three-contrib-log-link.bam')
    fasta = data_file('default-panel.fasta.gz')
    tmp_bam = str(tmp_path / 'three-contrib-log-link.bam')
    tmp_fasta = str(tmp_path / 'default-panel.fasta.gz')
    copyfile(bam, tmp_bam)
    copyfile(fasta, tmp_fasta)
    gt = microhapulator.type.type(tmp_bam, tmp_fasta)
    ac30 = gt.data['markers']['MHDBL000030']['allele_counts']
    ac197 = gt.data['markers']['MHDBL000197']['allele_counts']
    assert ac30 == {'A,A,T,C': 3, 'A,C,C,C': 2, 'A,C,C,G': 18, 'G,C,C,C': 1, 'G,C,C,G': 34}
    assert ac197 == {'A,A,T,T,T': 30, 'A,A,T,T,C': 39, 'A,A,T,A,T': 1, 'A,A,T,A,C': 1}


def test_type_cli_simple():
    with NamedTemporaryFile() as outfile:
        arglist = [
            'type', '--out', outfile.name, '--static', '5',
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

    gt = microhapulator.type.type(bam, fasta, static=10)
    assert gt.alleles('MHDBL000018') == set(['C,A,C,T,G', 'T,G,C,T,G'])
    assert gt.alleles('MHDBL000156') == set(['T,C,A,C', 'T,C,G,G'])

    gt = microhapulator.type.type(bam, fasta, static=4)
    assert gt.alleles('MHDBL000018') == set(['C,A,C,T,G', 'T,G,C,T,G', 'C,A,C,T,A', 'T,G,C,T,A'])
    assert gt.alleles('MHDBL000156') == set(['T,C,A,C', 'T,C,G,G'])


def test_type_no_var_offsets():
    bam = data_file('sandawe-dad.bam')
    fasta = data_file('sandawe-panel.fasta.gz')
    with pytest.raises(ValueError, match=r'variant offsets not annotated for target amplicon:'):
        profile = microhapulator.type.type(bam, fasta)
