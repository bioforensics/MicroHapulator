#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
from math import ceil
import microhapulator
import sys


def contrib(bamfile=None, refrfasta=None, pjson=None):
    if not pjson and (not bamfile or not refrfasta):
        message = 'must provide either JSON profile or BAM and refr FASTA'
        raise ValueError(message)
    if pjson:
        profile = microhapulator.profile.Profile(fromfile=pjson)
    else:
        profile = microhapulator.type.type(bamfile, refrfasta)
    num_alleles_per_marker = [len(profile.alleles(marker)) for marker in profile.markers()]
    max_num_alleles = max(num_alleles_per_marker)
    max_thresh = max_num_alleles - 1 if max_num_alleles % 2 == 0 else max_num_alleles
    max_loci = sum([1 for n in num_alleles_per_marker if n >= max_thresh])
    max_perc = round(max_loci / len(num_alleles_per_marker), 4)
    return ceil(max_num_alleles / 2), max_loci, max_perc


def main(args):
    ncontrib, nloci, ploci = contrib(bamfile=args.bam, refrfasta=args.refr, pjson=args.json)
    data = {
        'min_num_contrib': ncontrib,
        'num_loci_max_alleles': nloci,
        'perc_loci_max_alleles': ploci,
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
