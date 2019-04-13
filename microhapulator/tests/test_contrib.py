#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
from microhapulator.tests import data_file
import pytest


@pytest.mark.parametrize('gtjson,numcontrib', [
    ('single-contrib-1.json', 1),
    ('single-contrib-2.json', 1),
    ('single-contrib-3.json', 1),
    ('two-contrib-even.json', 2),
    ('three-contrib-even.json', 3),
    ('three-contrib-log.json', 3),
])
def test_contrib_json(gtjson, numcontrib):
    n, *data = microhapulator.contrib.contrib(gtjson=data_file(gtjson))
    assert n == numcontrib


def test_contrib_bam():
    bam = data_file('three-contrib-log.bam')
    refr = data_file('default-panel.fasta.gz')
    n, *data = microhapulator.contrib.contrib(bamfile=bam, refrfasta=refr)
    assert n == 3
