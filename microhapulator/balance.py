# -----------------------------------------------------------------------------
# Copyright (c) 2021, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

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
            readcount += mdata['num_discarded_reads']
        for haplotype, count in mdata['allele_counts'].items():
            readcount += count
        counts['Marker'].append(marker)
        counts['ReadCount'].append(readcount)
    data = pandas.DataFrame(counts)\
        .sort_values(['ReadCount'], ascending=False)\
        .reset_index(drop=True)
    return data


def balance(profile, include_discarded=True):
    data = count_and_sort(profile, include_discarded=include_discarded)
    with TemporaryDirectory() as tempdir:
        tfile = os.path.join(tempdir, "data.tsv")
        data.to_csv(tfile, index=False, header=False)
        subprocess.run(['termgraph', tfile])
    return data


def main(args):
    profile = microhapulator.profile.ObservedProfile(fromfile=args.input)
    data = balance(profile, include_discarded=args.discarded)
    if args.csv:
        data.to_csv(args.csv, index=False)
