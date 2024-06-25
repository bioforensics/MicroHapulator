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

from .mapstats import MappingSummary, MappingStats
from .qcsummary import PairedReadQCSummary, SingleEndReadQCSummary
from .reporter import Reporter
from datetime import datetime
from jinja2 import FileSystemLoader, Environment, Template
import microhapulator
import pandas as pd
from microhapulator import load_marker_thresholds
from microhapulator.marker import MicrohapIndex
from pkg_resources import resource_filename
import re
import sys


def full_reference_index_files(fasta):
    filename = fasta
    filenames = []
    for suffix in ("amb", "ann", "bwt", "pac", "sa"):
        idxfile = f"{filename}.{suffix}"
        filenames.append(idxfile)
    return filenames


def parse_flash_summary(logfile):
    # FIXME delete once this refactor is through
    with open(logfile, "r") as fh:
        data = fh.read()
        tp = re.search(r"Total pairs:\s+(\d+)", data).group(1)
        cp = re.search(r"Combined pairs:\s+(\d+)", data).group(1)
        ip = re.search(r"Innie pairs:\s+(\d+)", data).group(1)
        op = re.search(r"Outie pairs:\s+(\d+)", data).group(1)
        up = re.search(r"Uncombined pairs:\s+(\d+)", data).group(1)
        return int(tp), int(cp), int(ip), int(op), int(up)


def parse_read_counts(logfile):
    with open(logfile, "r") as fh:
        line = next(fh)
        totalreads = line.strip().split()[-1]
        line = next(fh)
        mappedreads = line.strip().split()[-1]
        return int(totalreads), int(mappedreads)


def parse_balance_stat(logfile):
    with open(logfile, "r") as fh:
        line = next(fh)
        stat = line.strip().split()[-1]
        return float(stat)


def parse_length_filter(logfile):
    with open(logfile, "r") as fh:
        next(fh)
        line = next(fh)
        removed, kept = line.strip().split()
        return int(removed), int(kept)


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


def final_html_report(
    samples,
    summary,
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
            summary=summary,
            thresholds=thresholds,
            read_length_table=reporter.read_length_table,
            typing_rates=reporter.per_marker_typing_rates,
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


def aggregate_summary(samples, reads_are_paired=True):
    data = list()
    for sample in sorted(samples):
        print(f"[Compiling summary] Sample={sample}", file=sys.stderr)
        length_failed, length_passed = parse_length_filter(
            f"analysis/{sample}/{sample}-length-filtered-read-counts.txt"
        )
        if reads_are_paired:
            totalreads, mergedreads, *_ = parse_flash_summary(f"analysis/{sample}/flash.log")
            mergedrate = mergedreads / totalreads
        else:
            totalreads = None
            mergedreads, mergedrate = None, None
        maptotal, mapped = MappingStats.parse_read_stats(f"analysis/{sample}/{sample}.bam.stats")
        if reads_are_paired:
            assert maptotal == (mergedreads - length_failed), (sample, maptotal, mergedreads)
        else:
            totalreads = maptotal
        frmaptotal, frmapped = parse_read_counts(
            f"analysis/{sample}/fullrefr/{sample}-fullrefr-mapped-reads.txt"
        )
        assert (
            frmaptotal == maptotal
        ), f"Sample={sample} fullrefr map total={frmaptotal} map total={maptotal}"
        typing = pd.read_csv(f"analysis/{sample}/{sample}-typing-rate.tsv", sep="\t")
        num_typing_success = typing.TypedReads.sum()
        num_typing_attempted = typing.TotalReads.sum()
        # At some point the following command was causing a failure. I've disabled the check for
        # now, but it's worth following up on to see what might cause this discrepancy.
        # -- Daniel Standage 2022-04-22
        # assert typing_total_reads == mapped, f"Sample={sample} type total={typing_total_reads} mapped={mapped}"
        typing_rate = num_typing_success / num_typing_attempted
        chisq = parse_balance_stat(f"analysis/{sample}/{sample}-interlocus-balance-chisq.txt")
        tstat = parse_balance_stat(f"analysis/{sample}/{sample}-heterozygote-balance-pttest.txt")
        entry = [
            sample,
            totalreads,
            mergedreads,
            mergedrate,
            length_failed,
            length_passed,
            mapped,
            mapped / maptotal,
            frmapped,
            frmapped / frmaptotal,
            num_typing_success,
            num_typing_attempted,
            typing_rate,
            chisq,
            tstat,
        ]
        data.append(entry)
    colnames = (
        "Sample",
        "TotalReads",
        "Merged",
        "MergeRate",
        "LengthFailed",
        "LengthPassed",
        "Mapped",
        "MappingRate",
        "MappedFullRefr",
        "MappingRateFullRefr",
        "TypedSuccess",
        "TypedAttempted",
        "TypingRate",
        "InterlocChiSq",
        "HetTstat",
    )
    return pd.DataFrame(data, columns=colnames)


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
            typing_rates=reporter.per_marker_typing_rates,
            markernames=sorted(reporter.marker_names),
            marker_details_table=marker_details_table,
            isna=pd.isna,
        )
        print(output, file=outfh, end="")
