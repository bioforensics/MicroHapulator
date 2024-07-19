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

from .details import MarkerDetails
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
    index = MicrohapIndex.from_files("marker-definitions.tsv", "marker-refr.fasta")
    marker_details = list(MarkerDetails.from_index(index))
    templatefile = resource_filename("microhapulator", "data/marker_details_template.html")
    with open(templatefile, "r") as infh, open("marker-detail-report.html", "w") as outfh:
        template = Template(infh.read())
        output = template.render(
            date=datetime.now().replace(microsecond=0).isoformat(),
            mhpl8rversion=microhapulator.__version__,
            mapping_rates=reporter.per_marker_mapping_rates,
            typing_summary=reporter.typing_summary,
            markernames=sorted(reporter.marker_names),
            marker_details=marker_details,
            isna=pd.isna,
        )
        print(output, file=outfh, end="")
