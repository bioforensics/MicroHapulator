# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

from .load import open, load_marker_frequencies
from .marker import MicrohapIndex
from .thresholds import ThresholdIndex
from . import profile
from . import api
from . import cli
