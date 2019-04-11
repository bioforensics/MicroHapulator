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
import pytest


def test_microhapulator_open():
    thefile = data_file('three-loci-refr.fasta')
    filehandle = microhapulator.open(thefile, 'r')
    filecontents = filehandle.read()
    assert len(filecontents.strip().split('\n')) == 6

    with pytest.raises(ValueError) as ve:
        filehandle = microhapulator.open(thefile, 'p')
    assert 'invalid mode "p"' in str(ve)
