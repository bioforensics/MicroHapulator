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
from contextlib import contextmanager
from gzip import open as gzopen
import os
from pkg_resources import resource_filename
import sys

# Internal modules
from microhapulator import profile

# Subcommands and command-line interface
from microhapulator import contain
from microhapulator import contrib
from microhapulator import dist
from microhapulator import diff
from microhapulator import mix
from microhapulator import prob
from microhapulator import seq
from microhapulator import sim
from microhapulator import type
from microhapulator import unite
from microhapulator import __main__
from microhapulator import cli


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


logstream = None
teelog = False


def package_file(path):
    return resource_filename('microhapulator', path)


@contextmanager
def open(filename, mode):
    if mode not in ('r', 'w'):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:
        filehandle = sys.stdin if mode == 'r' else sys.stdout
        yield filehandle
    else:
        openfunc = builtins.open
        if filename.endswith('.gz'):
            openfunc = gzopen
            mode += 't'
        with openfunc(filename, mode) as filehandle:
            yield filehandle


def plog(*args, **kwargs):
    """Print logging output."""
    global logstream
    global teelog
    if logstream is not None:
        print(*args, **kwargs, file=logstream)
    if logstream is None or teelog:
        print(*args, **kwargs, file=sys.stderr)
