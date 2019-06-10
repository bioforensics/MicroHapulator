#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
import microhapulator
from microhapulator.genotype import Genotype


def prob(popid, gt1, gt2=None):
    if gt2 is None:
        return gt1.rand_match_prob(popid)
    else:
        return gt1.rmp_lr_test(gt2, popid)


def main(args):
    gt1 = Genotype(fromfile=args.genotype1)
    gt2 = Genotype(fromfile=args.genotype2) if args.genotype2 else None
    result = prob(args.population, gt1, gt2)
    key = 'random_match_probability' if gt2 is None else 'rmp_likelihood_ratio'
    data = {
        key: '{:.3E}'.format(result),
    }
    with microhapulator.open(args.out, 'w') as fh:
        json.dump(data, fh, indent=4)
