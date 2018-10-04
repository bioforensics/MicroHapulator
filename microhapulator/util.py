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
import pyfaidx


def data_file(path):
    pathparts = path.split('/')
    relpath = os.path.join('data', *pathparts)
    return resource_filename('microhapulator', relpath)


bogus_loci = [
    {
        'ID': 'FakeLocus1',
        'Name': 'mh18DS-42',
        'Chrom': 2,
        'Start': 1033,
        'End': 1234,
        'Variants': 'A,B,C',
    },
]

bogus_index = pyfaidx.Fasta(data_file('bogus-refr.fa'))
