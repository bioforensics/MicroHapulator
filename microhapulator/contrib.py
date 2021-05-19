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
from microhapulator.profile import Profile
import sys


def load_profile(bamfile=None, refrfasta=None, json=None, **kwargs):
    if not json and (not bamfile or not refrfasta):
        message = 'must provide either JSON profile or BAM and refr FASTA'
        raise ValueError(message)
    if json:
        profile = microhapulator.profile.Profile(fromfile=json)
    else:
        profile = microhapulator.type.type(bamfile, refrfasta, **kwargs)
    return profile


def contrib(profile):
    num_alleles_per_marker = [len(profile.alleles(marker)) for marker in profile.markers()]
    max_num_alleles = max(num_alleles_per_marker)
    max_thresh = max_num_alleles - 1 if max_num_alleles % 2 == 0 else max_num_alleles
    max_loci = sum([1 for n in num_alleles_per_marker if n >= max_thresh])
    max_perc = round(max_loci / len(num_alleles_per_marker), 4)
    return ceil(max_num_alleles / 2), max_loci, max_perc


def main(args):
    profile = load_profile(
        bamfile=args.bam, refrfasta=args.refr, json=args.json, static=args.static,
        dynamic=args.dynamic
    )
    ncontrib, nloci, ploci = contrib(profile)
    data = {
        'min_num_contrib': ncontrib,
        'num_loci_max_alleles': nloci,
        'perc_loci_max_alleles': ploci,
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
