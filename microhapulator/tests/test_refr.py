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
    seqiter = microhapulator.refr.get_seqs(['MHDBL000185'], seqindex)
    seqs = list(seqiter)
    assert len(seqs) == 1
    assert seqs[0][0] == 'MHDBL000185 GRCh38:chr5:38881266-38881617 variants=91:169:234:259'
    assert seqs[0][1].startswith('CTGGCACAGTGAGCACCTTCTGTCTCTGATCTGTT')
    assert seqs[0][1].endswith('TTGCTTTTAGGGGAATTACAGCACCACTGTGAAGT')


def test_refr_cli_simple():
    with NamedTemporaryFile() as outfile:
        arglist = [
            'refr', '--min-length', '250', '--out', outfile.name, 'hg38.fasta',
            'mh06KK-008', 'SI664615C', 'MHDBL000031'
        ]
        microhapulator.__main__.main(arglist)
        testoutfile = data_file('three-loci-refr.fasta')
        assert filecmp.cmp(outfile.name, testoutfile)
