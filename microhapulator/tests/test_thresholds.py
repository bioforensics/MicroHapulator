# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from microhapulator import ThresholdIndex
from microhapulator.tests import data_file
import pytest


def test_load_marker_filters():
    configfile = data_file("filters.csv")
    thresholds = ThresholdIndex.load(configfile=configfile, global_static=10, global_dynamic=0.015)
    assert thresholds.get("mh01XYZ-1") == (5, pytest.approx(0.001))
    assert thresholds.get("mh01XYZ-2") == (10, pytest.approx(0.015))
    assert thresholds.get("mh01XYZ-3") == (10, pytest.approx(0.015))


def test_load_marker_filters_default():
    thresholds = ThresholdIndex.load()
    for marker in ("mh04ZZZ-1", "mh06ZZZ-2", "mh19ZZZ-3", "mh22ZZZ-4"):
        static, dynamic = thresholds.get(marker)
        assert static == 5
        assert dynamic == pytest.approx(0.02)
