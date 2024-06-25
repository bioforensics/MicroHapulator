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
from datetime import datetime
from jinja2 import FileSystemLoader, Environment, Template
import microhapulator
import pandas as pd
from microhapulator import load_marker_thresholds
from microhapulator.marker import MicrohapIndex
from pkg_resources import resource_filename


def final_html_report(
    samples,
    reads_are_paired=True,
    thresh_static=10,
    thresh_dynamic=0.02,
    thresh_file=None,
    ambiguous_read_threshold=0.2,
    length_threshold=50,
):
    reporter = Reporter(samples, reads_are_paired=reads_are_paired)
    thresholds = load_marker_thresholds(
        reporter.marker_names, thresh_file, thresh_static, thresh_dynamic
    )
    thresholds.fillna(0, inplace=True)
    template_loader = FileSystemLoader(resource_filename("microhapulator", "data"))
    env = Environment(loader=template_loader)
    if reads_are_paired:
        template_file = "paired.html"
    else:
        template_file = "single.html"
    template = env.get_template(template_file)
    with open("report.html", "w") as outfh:
        output = template.render(
            date=datetime.now().replace(microsecond=0).isoformat(),
            mhpl8rversion=microhapulator.__version__,
            samples=samples,
            thresholds=thresholds,
            read_length_table=reporter.read_length_table,
            typing_summary=reporter.typing_summary,
            mapping_rates=reporter.per_marker_mapping_rates,
            markernames=reporter.marker_names,
            qc=reporter.qc_summary,
            mapping_summary=reporter.mapping_summary,
            repetitive_reads_by_marker=reporter.mapping_summary.repetitive_reads_by_marker(),
            reads_are_paired=reads_are_paired,
            ambiguous_read_threshold=ambiguous_read_threshold,
            read_length_threshold=length_threshold,
        )
        print(output, file=outfh, end="")


def marker_detail_report(samples, reads_are_paired=True):
    reporter = Reporter(samples, reads_are_paired=reads_are_paired)
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
