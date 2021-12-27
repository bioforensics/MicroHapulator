# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import microhapulator
from microhapulator.profile import SimulatedProfile


def main(args):
    profiles = [SimulatedProfile(pfile) for pfile in args.profiles]
    combined = SimulatedProfile.merge(profiles)
    with microhapulator.open(args.out, "w") as fh:
        combined.dump(fh)
