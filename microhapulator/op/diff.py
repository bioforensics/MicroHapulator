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


def diff(prof1, prof2):
    allmarkers = set(prof1.markers()).union(prof2.markers())
    for marker in sorted(allmarkers):
        allele1 = prof1.alleles(marker)
        allele2 = prof2.alleles(marker)
        diff1 = allele1 - allele2
        diff2 = allele2 - allele1
        if len(diff1) > 0 or len(diff2) > 0:
            yield marker, diff1, diff2
