#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import os
from pkg_resources import resource_filename


def data_file(path):
    pathparts = path.split('/')
    relpath = os.path.join('data', *pathparts)
    return resource_filename('microhapulator', relpath)
