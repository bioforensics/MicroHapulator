#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
import microhapdb
import microhapulator
from microhapulator.profile import Profile


def prob(popid, prof1, prof2=None, erate=0.001):
    freqs = microhapdb.frequencies[microhapdb.frequencies.Population == popid]
    if prof2 is None:
        return prof1.rand_match_prob(freqs)
    else:
        return prof1.rmp_lr_test(prof2, freqs, erate=erate)


def main(args):
    prof1 = Profile(fromfile=args.profile1)
    prof2 = Profile(fromfile=args.profile2) if args.profile2 else None
    popids = microhapdb.population.standardize_ids([args.population])
    if len(popids) != 1:
        message = 'issue with population "{:s}"; invalid or not unique'.format(args.population)
        raise ValueError(message)
    popid = popids.iloc[0]
    result = prob(popid, prof1, prof2=prof2, erate=args.erate)
    key = "random_match_probability" if prof2 is None else "likelihood_ratio"
    data = {
        key: "{:.3E}".format(result),
    }
    with microhapulator.open(args.out, "w") as fh:
        json.dump(data, fh, indent=4)
