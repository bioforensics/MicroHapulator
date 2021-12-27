# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------

import microhapulator
from microhapulator.tests import data_file
import pytest


def test_microhapulator_open():
    thefile = data_file("refr/three-loci-refr.fasta")
    with microhapulator.open(thefile, "r") as filehandle:
        filecontents = filehandle.read()
        assert len(filecontents.strip().split("\n")) == 6
    with pytest.raises(ValueError, match=r'invalid mode "p"') as ve:
        with microhapulator.open(thefile, "p") as filehandle:
            pass
