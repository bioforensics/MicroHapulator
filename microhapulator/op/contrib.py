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

from math import ceil


def contrib(profile):
    num_alleles_per_marker = [len(profile.alleles(marker)) for marker in profile.markers()]
    max_num_alleles = max(num_alleles_per_marker)
    max_thresh = max_num_alleles - 1 if max_num_alleles % 2 == 0 else max_num_alleles
    max_loci = sum([1 for n in num_alleles_per_marker if n >= max_thresh])
    max_perc = round(max_loci / len(num_alleles_per_marker), 4)
    return ceil(max_num_alleles / 2), max_loci, max_perc
