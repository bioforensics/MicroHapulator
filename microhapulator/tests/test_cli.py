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
    with microhapulator.open(thefile, 'r') as filehandle:
        filecontents = filehandle.read()
        assert len(filecontents.strip().split('\n')) == 6

    with pytest.raises(ValueError, match=r'invalid mode "p"') as ve:
        with microhapulator.open(thefile, 'p') as filehandle:
            pass
