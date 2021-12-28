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

from collections import defaultdict
import microhapulator
import os.path
import pysam
import re


def check_index(bamfile):
    index1 = bamfile + ".bai"
    index2 = re.sub(r".bam$", ".bai", bamfile)
    for testfile in (index1, index2):
        if os.path.isfile(testfile):
            break
    else:
        pysam.index(bamfile)


def tally_haplotypes(bam, offsets, minbasequal=10, max_depth=1e6):
    totaldiscarded = 0
    for locusid in sorted(offsets):
        discarded = 0
        haplotypes = defaultdict(int)
        ht = defaultdict(dict)
        varloc = set(offsets[locusid])
        cov_pos = list()
        for column in bam.pileup(locusid, min_base_quality=minbasequal, max_depth=max_depth):
            cov_pos.append(column.n)
            if column.pos not in varloc:
                continue
            for record in column.pileups:
                aligned_base = None
                if record.is_del or record.is_refskip:
                    continue
                aligned_base = record.alignment.query_sequence[record.query_position]
                ht[record.alignment.query_name][column.pos] = aligned_base
        for readname, htdict in ht.items():
            htlist = [htdict[pos] for pos in sorted(htdict)]
            if len(htlist) < len(varloc):
                discarded += 1
                continue
            htstr = ",".join(htlist)
            haplotypes[htstr] += 1
        yield locusid, cov_pos, haplotypes, discarded
        totaldiscarded += discarded
    microhapulator.plog(
        "[MicroHapulator::type] discarded",
        totaldiscarded,
        "reads with gaps or missing data at positions of interest",
    )


def type(
    bamfile, markertsv, minbasequal=10, ecthreshold=0.25, static=None, dynamic=None, max_depth=1e6
):
    check_index(bamfile)
    bam = pysam.AlignmentFile(bamfile, "rb")
    markers = microhapulator.load_marker_definitions(markertsv)
    offsets = defaultdict(list)
    for n, row in markers.iterrows():
        offsets[row.Marker].append(row.Offset)
    microhapulator.cross_check_marker_ids(
        bam.references, offsets.keys(), "read alignments", "marker definitions"
    )
    genotyper = tally_haplotypes(bam, offsets, minbasequal=minbasequal, max_depth=max_depth)
    profile = microhapulator.profile.ObservedProfile()
    for locusid, cov_by_pos, htcounts, ndiscarded in genotyper:
        profile.record_coverage(locusid, cov_by_pos, ndiscarded=ndiscarded)
        for allele, count in htcounts.items():
            profile.record_allele(locusid, allele, count)
    profile.infer(ecthreshold=ecthreshold, static=static, dynamic=dynamic)
    return profile
