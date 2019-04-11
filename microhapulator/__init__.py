#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from sys import stderr

from microhapulator import locus
from microhapulator import population

from microhapulator import sim
from microhapulator import cli
from microhapulator.util import data_file, bogus_loci, bogus_index
from microhapulator.context import LocusContext
from microhapulator.genotype import Genotype


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


logstream = None
teelog = False


def plog(*args, **kwargs):
    """Print logging output."""
    if logstream is not None:
        print(*args, **kwargs, file=logstream)
    if logstream is None or teelog:
        print(*args, **kwargs, file=stderr)
