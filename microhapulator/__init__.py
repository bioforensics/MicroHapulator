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

import builtins
from contextlib import contextmanager
from gzip import open as gzopen
import pandas as pd
import sys
from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


@contextmanager
def open(filename, mode, enter=True):
    if mode not in ("r", "w"):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ["-", None]:
        filehandle = sys.stdin if mode == "r" else sys.stdout
        yield filehandle
    else:
        openfunc = builtins.open
        if filename.endswith(".gz"):
            openfunc = gzopen
            mode += "t"
        if enter:
            with openfunc(filename, mode) as filehandle:
                yield filehandle
        else:
            yield openfunc(filename, mode)


def load_marker_frequencies(tsvfile):
    frequencies = pd.read_csv(tsvfile, sep="\t")
    missing = set(["Marker", "Haplotype", "Frequency"]) - set(frequencies.columns)
    if len(missing) > 0:
        message = "column(s) missing from marker frequency file: " + ", ".join(sorted(missing))
        raise ValueError(message)
    return frequencies


def load_marker_thresholds(markernames, configfile=None, global_static=5, global_dynamic=0.02):
    default = (global_static, global_dynamic)
    thresholds = pd.DataFrame([default], index=sorted(markernames), columns=["Static", "Dynamic"])
    if configfile:
        config = pd.read_csv(configfile, sep=None, engine="python")
        missing = set(["Marker", "Static", "Dynamic"]) - set(config.columns)
        if len(missing) > 0:
            missingstr = ",".join(sorted(missing))
            raise ValueError(f"filter config file missing column(s): {missingstr}")
        if len(config.Marker) != len(config.Marker.unique()):
            raise ValueError("filter config file contains duplicate entries for some markers")
        config = config.set_index("Marker", drop=True)
        for marker, row in config.iterrows():
            thresholds.loc[marker, "Static"] = row.Static
            thresholds.loc[marker, "Dynamic"] = row.Dynamic
    return thresholds.astype({"Static": "Int64"})


from .marker import MicrohapIndex
from . import profile
from . import api
from . import cli
from . import pipeaux
