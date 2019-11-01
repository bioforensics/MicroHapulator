#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
import numpy.random
import microhapulator
from microhapulator.profile import Profile


def main(args):
    if args.seed:
        numpy.random.seed(args.seed)
    profile = Profile.unite(
        Profile(fromfile=args.mom),
        Profile(fromfile=args.dad),
    )
    with microhapulator.open(args.out, 'w') as fh:
        profile.dump(fh)
