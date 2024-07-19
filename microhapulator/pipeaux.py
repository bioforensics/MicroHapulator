# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from .reporter import Reporter
from .thresholds import ThresholdIndex
from datetime import datetime
from jinja2 import Template
import microhapulator
import pandas as pd
from microhapulator.marker import MicrohapIndex
from pkg_resources import resource_filename


def marker_detail_report(samples, reads_are_paired=True):
    reporter = Reporter(samples, ThresholdIndex(), reads_are_paired=reads_are_paired)
    marker_details_table = marker_details()
    templatefile = resource_filename("microhapulator", "data/marker_details_template.html")
    with open(templatefile, "r") as infh, open("marker-detail-report.html", "w") as outfh:
        template = Template(infh.read())
        output = template.render(
            date=datetime.now().replace(microsecond=0).isoformat(),
            mhpl8rversion=microhapulator.__version__,
            mapping_rates=reporter.per_marker_mapping_rates,
            typing_summary=reporter.typing_summary,
            markernames=sorted(reporter.marker_names),
            marker_details_table=marker_details_table,
            isna=pd.isna,
        )
        print(output, file=outfh, end="")


def marker_details():
    index = MicrohapIndex.from_files("marker-definitions.tsv", "marker-refr.fasta")
    all_marker_details = list()
    for locus, marker in index:
        marker_offsets = ", ".join([str(o) for o in sorted(marker.offsets_locus)])
        chrom = marker.chrom
        offsets38 = ", ".join([str(o) for o in sorted(marker.offsets_chrom)])
        seq = locus.sequence.strip().upper()
        gc_content = round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)
        sample_details = [
            marker.id,
            len(seq),
            gc_content,
            marker_offsets,
            seq,
            chrom,
            offsets38,
        ]
        all_marker_details.append(sample_details)
    col_names = ["Marker", "Length", "GC", "Offsets", "Sequence", "Chrom", "Hg38Offset"]
    marker_details_table = pd.DataFrame(all_marker_details, columns=col_names).set_index("Marker")
    return marker_details_table
