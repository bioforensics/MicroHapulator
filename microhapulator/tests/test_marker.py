# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import microhapulator
from microhapulator.marker import DefinitionIndex
from microhapulator.tests import data_file
import pytest


def test_index_basic():
    defn = data_file("def/2x2.tsv")
    index = DefinitionIndex.from_csv(defn, strict=True)
    assert index.loci == ["mh21FHL-002", "mh21KK-320"]
    assert index["mh21KK-320"]["mh21KK-320.v1"] == [10, 80, 169, 195]
