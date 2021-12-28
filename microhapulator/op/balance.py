# -------------------------------------------------------------------------------------------------
# Copyright (c) 2021, DHS.
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
import os
import pandas
from tempfile import TemporaryDirectory
import subprocess


def count_and_sort(profile, include_discarded=True):
    counts = dict(
        Marker=list(),
        ReadCount=list(),
    )
    for marker, mdata in profile.data["markers"].items():
        readcount = 0
        if include_discarded:
            readcount += mdata["num_discarded_reads"]
        for haplotype, count in mdata["allele_counts"].items():
            readcount += count
        counts["Marker"].append(marker)
        counts["ReadCount"].append(readcount)
    data = (
        pandas.DataFrame(counts).sort_values(["ReadCount"], ascending=False).reset_index(drop=True)
    )
    return data


def balance(profile, include_discarded=True):
    data = count_and_sort(profile, include_discarded=include_discarded)
    with TemporaryDirectory() as tempdir:
        tfile = os.path.join(tempdir, "data.tsv")
        data.to_csv(tfile, index=False, header=False)
        subprocess.run(["termgraph", tfile])
    return data
