# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
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
import microhapulator.api as mhapi
from microhapulator.profile import TypingResult
from microhapulator.tests import data_file
import pandas as pd
import pytest


@pytest.mark.parametrize(
    "counts,expfile",
    [
        (True, "profile-efm.csv"),
        (False, "profile-lrmix.csv"),
    ],
)
def test_convert_counts(tmp_path, counts, expfile):
    csvfile = str(tmp_path / "out.csv")
    result = TypingResult(fromfile=data_file("prof/deep-filt-clean.json"))
    result.dump_csv(csvfile, "MySample", counts=counts)
    observed = pd.read_csv(csvfile)
    expected = pd.read_csv(data_file(expfile))
    assert observed.equals(expected)
