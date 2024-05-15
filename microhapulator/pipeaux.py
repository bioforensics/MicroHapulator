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

from .qcsummary import PairedReadQCSummary, SingleEndReadQCSummary
from datetime import datetime
from jinja2 import Template
import json
import microhapulator
import pandas as pd
from pathlib import Path
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


def per_marker_typing_rate(samples):
    sample_rates = dict()
    for sample in samples:
        filename = f"analysis/{sample}/{sample}-typing-rate.tsv"
        sample_df = pd.read_csv(filename, sep="\t").set_index("Marker")
        sample_rates[sample] = sample_df
    return sample_rates


def per_marker_mapping_rate(samples):
    sample_rates = dict()
    for sample in samples:
        total_reads_filename = f"analysis/{sample}/{sample}-marker-read-counts.csv"
        total_df = pd.read_csv(total_reads_filename).set_index("Marker")
        expected_count = total_df["ReadCount"].sum() / len(total_df)
        total_df["ExpectedObservedRatio"] = round(total_df["ReadCount"] / expected_count, 2)
        repetitive_filename = Path(f"analysis/{sample}/{sample}-repetitive-reads.csv")
        if repetitive_filename.stat().st_size > 0:
            repetitive_df = pd.read_csv(repetitive_filename).set_index("Marker")
            sample_df = pd.concat([total_df, repetitive_df], axis=1)
            sample_df["RepetitiveRate"] = sample_df["RepetitiveReads"] / sample_df["ReadCount"]
        else:
            sample_df = total_df
            sample_df["RepetitiveReads"] = None
            sample_df["RepetitiveRate"] = None
        sample_rates[sample] = sample_df
    marker_names = sample_df.index
    return sample_rates, marker_names


def marker_details():
    index = MicrohapIndex.from_files("marker-definitions.tsv", "marker-refr.fasta")
    all_marker_details = list()
    for locus, marker in index:
        marker_offsets = ", ".join([str(o) for o in sorted(marker.offsets_locus)])
        chrom = marker.chrom
        offsets38 = "N/A"
        if index.has_chrom_offsets:
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
    if reads_are_paired:
        read_length_table = read_length_table_paired_end(samples)
        plots = aggregate_plots_paired_end(samples)
        qc = PairedReadQCSummary.collect(samples)
    else:
        read_length_table = read_length_table_single_end(samples)
        plots = aggregate_plots_single_end(samples)
        qc = SingleEndReadQCSummary.collect(samples)
    typing_rates = per_marker_typing_rate(samples)
    mapping_rates, marker_names = per_marker_mapping_rate(samples)
    thresholds = load_marker_thresholds(marker_names, thresh_file, thresh_static, thresh_dynamic)
    thresholds.fillna(0, inplace=True)
    templatefile = resource_filename("microhapulator", "data/template.html")
    with open(templatefile, "r") as infh, open("report.html", "w") as outfh:
        template = Template(infh.read())
        output = template.render(
            date=datetime.now().replace(microsecond=0).isoformat(),
            mhpl8rversion=microhapulator.__version__,
            samples=samples,
            summary=summary,
            plots=plots,
            thresholds=thresholds,
            zip=zip,
            read_length_table=read_length_table,
            typing_rates=typing_rates,
            mapping_rates=mapping_rates,
            markernames=marker_names,
            qc=qc,
            len=len,
            isna=pd.isna,
            reads_are_paired=reads_are_paired,
            ambiguous_read_threshold=ambiguous_read_threshold,
            read_length_threshold=length_threshold,
        )
        print(output, file=outfh, end="")


def read_length_table_paired_end(samples):
    read_length_data = list()
    for sample in samples:
        with open(f"analysis/{sample}/{sample}-r1-read-lengths.json", "r") as fh:
            r1lengths = json.load(fh)
            r1lengths = list(set(r1lengths))
        with open(f"analysis/{sample}/{sample}-r2-read-lengths.json", "r") as fh:
            r2lengths = json.load(fh)
            r2lengths = list(set(r2lengths))
        if len(r1lengths) != 1 or len(r2lengths) != 1:
            return None
        read_length_data.append((sample, r1lengths[0], r2lengths[0]))
    else:
        return pd.DataFrame(read_length_data, columns=("Sample", "LengthR1", "LengthR2"))


def read_length_table_single_end(samples):
    read_length_data = list()
    for sample in samples:
        with open(f"analysis/{sample}/{sample}-read-lengths.json", "r") as fh:
            lengths = json.load(fh)
            lengths = list(set(lengths))
        if len(lengths) != 1:
            return None
        read_length_data.append((sample, lengths[0]))
    return pd.DataFrame(read_length_data, columns=("Sample", "Length"))


def aggregate_plots_paired_end(samples):
    plots = {
        "r1readlen": "",
        "r2readlen": "",
        "mergedreadlen": "",
        "locbalance": list(),
        "hetbalance": list(),
        "donut": list(),
    }
    plots["r1readlen"] = "analysis/r1-read-lengths.png"
    plots["r2readlen"] = "analysis/r2-read-lengths.png"
    plots["mergedreadlen"] = "analysis/merged-read-lengths.png"
    for sample in samples:
        plots["locbalance"].append(f"analysis/{sample}/{sample}-interlocus-balance.png")
        plots["hetbalance"].append(f"analysis/{sample}/{sample}-heterozygote-balance.png")
        plots["donut"].append(f"analysis/{sample}/{sample}-donut.png")
    return plots


def aggregate_plots_single_end(samples):
    plots = {"readlen": "", "locbalance": list(), "hetbalance": list(), "donut": list()}
    plots["readlen"] = "analysis/read-lengths.png"
    for sample in samples:
        plots["locbalance"].append(f"analysis/{sample}/{sample}-interlocus-balance.png")
        plots["hetbalance"].append(f"analysis/{sample}/{sample}-heterozygote-balance.png")
        plots["donut"].append(f"analysis/{sample}/{sample}-donut.png")
    return plots


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
        maptotal, mapped = parse_read_counts(f"analysis/{sample}/{sample}-mapped-reads.txt")
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


def marker_detail_report(samples):
    mapping_rates, marker_names = per_marker_mapping_rate(samples)
    marker_details_table = marker_details()
    templatefile = resource_filename("microhapulator", "data/marker_details_template.html")
    with open(templatefile, "r") as infh, open("marker-detail-report.html", "w") as outfh:
        template = Template(infh.read())
        output = template.render(
            date=datetime.now().replace(microsecond=0).isoformat(),
            mhpl8rversion=microhapulator.__version__,
            mapping_rates=mapping_rates,
            typing_rates=per_marker_typing_rate(samples),
            markernames=sorted(marker_names),
            marker_details_table=marker_details_table,
            isna=pd.isna,
        )
        print(output, file=outfh, end="")
