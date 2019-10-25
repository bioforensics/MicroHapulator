#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapulator
from microhapulator.profile import SimulatedProfile


def main(args):
    genotypes = [SimulatedGenotype(gtfile) for gtfile in args.gtfiles]
    combined = SimulatedGenotype.merge(genotypes)
    with microhapulator.open(args.out, 'w') as fh:
        combined.dump(fh)
