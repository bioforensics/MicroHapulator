#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import filecmp
import microhapulator
from microhapulator.tests import data_file
import pyfaidx
import pytest
from tempfile import NamedTemporaryFile


def test_refr_simple():
    seqindex = pyfaidx.Fasta('hg38.fasta')
    seqiter = microhapulator.refr.get_seqs(['MHDBL000146'], seqindex)
    seqs = list(seqiter)
    assert len(seqs) == 1
    assert seqs[0][0] == 'MHDBL000146 GRCh38:chr22:44857743-44858094'
    assert seqs[0][1].startswith('CCAGCCTCCCACATGCAAGCTGTGTGACCTCGGGC')
    assert seqs[0][1].endswith('TTTCAAGGGAATTCCTTGTCCATTCAAATGACTGA')


def test_refr_cli_simple():
    with NamedTemporaryFile() as outfile:
        arglist = [
            'refr', '--min-length', '250', '--out', outfile.name, 'hg38.fasta',
            'mh06KK-008', 'SI664615C', 'MHDBL000031'
        ]
        microhapulator.__main__.main(arglist)
        testoutfile = data_file('three-loci-refr.fasta')
        assert filecmp.cmp(outfile.name, testoutfile)
