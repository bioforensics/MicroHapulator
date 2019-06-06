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


def contrib(bamfile=None, refrfasta=None, gtjson=None):
    if not gtjson and (not bamfile or not refrfasta):
        message = 'must provide either genotype JSON or BAM and refr FASTA'
        raise ValueError(message)
    if gtjson:
        gt = microhapulator.genotype.Genotype(fromfile=gtjson)
    else:
        gt = microhapulator.type.type(bamfile, refrfasta)
    num_alleles_per_locus = [len(gt.alleles(locusid)) for locusid in gt.loci()]
    max_num_alleles = max(num_alleles_per_locus)
    max_thresh = max_num_alleles - 1 if max_num_alleles % 2 == 0 else max_num_alleles
    max_loci = sum([1 for n in num_alleles_per_locus if n >= max_thresh])
    max_perc = round(max_loci / len(num_alleles_per_locus), 4)
    return ceil(max_num_alleles / 2), max_loci, max_perc


def main(args):
    ncontrib, nloci, ploci = contrib(bamfile=args.bam, refrfasta=args.refr, gtjson=args.json)
    data = {
        'min_num_contrib': ncontrib,
        'num_loci_max_alleles': nloci,
        'perc_loci_max_alleles': ploci,
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
