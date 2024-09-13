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


@contextmanager
def open(filename, mode):
    if mode not in ("r", "w"):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ["-", None]:
        filehandle = sys.stdin if mode == "r" else sys.stdout
        yield filehandle
    else:
        openfunc = builtins.open
        if str(filename).endswith(".gz"):
            openfunc = gzopen
            mode += "t"
        with openfunc(filename, mode) as filehandle:
            yield filehandle


def load_marker_frequencies(tsvfile):
    frequencies = pd.read_csv(tsvfile, sep="\t")
    missing = set(["Marker", "Haplotype", "Frequency"]) - set(frequencies.columns)
    if len(missing) > 0:
        message = "column(s) missing from marker frequency file: " + ", ".join(sorted(missing))
        raise ValueError(message)
    return frequencies
