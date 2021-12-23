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
from microhapulator.profile import Profile


def prob(frequencies, prof1, prof2=None, erate=0.001):
    if prof2 is None:
        return prof1.rand_match_prob(frequencies)
    else:
        return prof1.rmp_lr_test(prof2, frequencies, erate=erate)


def main(args):
    prof1 = Profile(fromfile=args.profile1)
    prof2 = Profile(fromfile=args.profile2) if args.profile2 else None
    frequencies = microhapulator.load_marker_frequencies(args.freq)
    result = prob(frequencies, prof1, prof2=prof2, erate=args.erate)
    key = "random_match_probability" if prof2 is None else "likelihood_ratio"
    data = {
        key: "{:.3E}".format(result),
    }
    with microhapulator.open(args.out, "w") as fh:
        json.dump(data, fh, indent=4)
