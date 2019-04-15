#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

# Core libraries
import builtins
from gzip import open as gzopen
from sys import stdin, stdout, stderr

# Internal modules
from microhapulator import locus
from microhapulator import population

# Subcommands and command-line interface
from microhapulator import contrib
from microhapulator import dist
from microhapulator import refr
from microhapulator import sim
from microhapulator import type
from microhapulator import __main__
from microhapulator import cli


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


logstream = None
teelog = False


def open(filename, mode):
    if mode not in ('r', 'w'):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:
        filehandle = stdin if mode == 'r' else stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


def plog(*args, **kwargs):
    """Print logging output."""
    if logstream is not None:
        print(*args, **kwargs, file=logstream)
    if logstream is None or teelog:
        print(*args, **kwargs, file=stderr)
