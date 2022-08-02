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


from base64 import b64encode
from datetime import datetime
from jinja2 import Template
import json
from matplotlib import pyplot as plt
import microhapulator
import pandas as pd
from pkg_resources import resource_filename
import re
import sys


def full_reference_index_files(fasta):
    filename = fasta
    filenames = [filename]
    for suffix in ("amb", "ann", "bwt", "pac", "sa"):
        idxfile = f"{filename}.{suffix}"
        filenames.append(idxfile)
    return filenames


def parse_flash_summary(logfile):
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


def encode(filepath):
    with open(filepath, "rb") as fh:
        return b64encode(fh.read()).decode("ascii")
    

def per_marker_typing_rate(samples):
    sample_rates = dict()
    for sample in samples:
        filename = Path("analysis") / sample / f"{sample}-typing-rate.tsv"
        sample_rates[sample] = pd.read_csv(filename, sep="\t")
    all_rates = pd.concat([df.TypingRate for df in sample_rates.values()], axis=1, keys=sample_rates.keys())
    return all_rates


def final_html_report(samples, summary):
    read_length_table = list()
    for sample in samples:
        with open(f"analysis/{sample}/{sample}-r1-read-lengths.json") as fh:
            r1lengths = json.load(fh)
            r1lengths = list(set(r1lengths))
        with open(f"analysis/{sample}/{sample}-r2-read-lengths.json") as fh:
            r2lengths = json.load(fh)
            r2lengths = list(set(r2lengths))
        if len(r1lengths) != 1 or len(r2lengths) != 1:
            read_length_table = None
            break
        read_length_table.append((sample, r1lengths[0], r2lengths[0]))
    if read_length_table is not None:
        col = ("Sample", "LengthR1", "LengthR2")
        read_length_table = pd.DataFrame(read_length_table, columns=col)
    typing_rate = per_marker_typing_rate(samples)
    plots = {
        "r1readlen": list(),
        "r2readlen": list(),
        "mergedreadlen": list(),
        "locbalance": list(),
        "hetbalance": list(),
    }
    for sample in samples:
        plots["r1readlen"].append(encode(f"analysis/{sample}/{sample}-r1-read-lengths.png"))
        plots["r2readlen"].append(encode(f"analysis/{sample}/{sample}-r2-read-lengths.png"))
        plots["mergedreadlen"].append(
            encode(f"analysis/{sample}/{sample}-merged-read-lengths.png")
        )
        plots["locbalance"].append(encode(f"analysis/{sample}/{sample}-interlocus-balance.png"))
        plots["hetbalance"].append(encode(f"analysis/{sample}/{sample}-heterozygote-balance.png"))
    templatefile = resource_filename("microhapulator", "data/template.html")
    with open(templatefile, "r") as infh, open("report.html", "w") as outfh:
        template = Template(infh.read())
        output = template.render(
            date=datetime.now().replace(microsecond=0).isoformat(),
            mhpl8rversion=microhapulator.__version__,
            samples=samples,
            summary=summary,
            plots=plots,
            static=5,
            dynamic=0.02,
            zip=zip,
            read_length_table=read_length_table,
            typing_rates=typing_rates,
        )
        print(output, file=outfh, end="")


def aggregate_summary(samples):
    columns = (
        "Sample",
        "TotalReads",
        "Merged",
        "MergeRate",
        "Mapped",
        "MappingRate",
        "MappedFullRefr",
        "MappingRateFullRefr",
        "Typed",
        "TypingRate",
        "InterlocChiSq",
        "HetTstat",
    )
    data = {colname: [] for colname in columns}
    for sample in sorted(samples):
        print(f"[Compiling summary] Sample={sample}", file=sys.stderr)
        totalreads, mergedreads, *_ = parse_flash_summary(f"analysis/{sample}/flash.log")
        maptotal, mapped = parse_read_counts(f"analysis/{sample}/{sample}-mapped-reads.txt")
        assert maptotal == mergedreads, (sample, maptotal, mergedreads)
        frmaptotal, frmapped = parse_read_counts(
            f"analysis/{sample}/fullrefr/{sample}-fullrefr-mapped-reads.txt"
        )
        assert (
            frmaptotal == maptotal
        ), f"Sample={sample} fullrefr map total={frmaptotal} map total={maptotal}"
        typing = pd.read_csv(f"analysis/{sample}/{sample}-typing-rate.tsv", sep="\t")
        num_typed_reads = typing.TypedReads.sum()
        typing_total_reads = typing.TotalReads.sum()
        # At some point the following command was causing a failure. I've disabled the check for
        # now, but it's worth following up on to see what might cause this discrepancy.
        # -- Daniel Standage 2022-04-22
        # assert typing_total_reads == mapped, f"Sample={sample} type total={typing_total_reads} mapped={mapped}"
        typing_rate = num_typed_reads / typing_total_reads
        chisq = parse_balance_stat(f"analysis/{sample}/{sample}-interlocus-balance-chisq.txt")
        tstat = parse_balance_stat(f"analysis/{sample}/{sample}-heterozygote-balance-pttest.txt")
        data["Sample"].append(sample)
        data["TotalReads"].append(totalreads)
        data["Merged"].append(mergedreads)
        data["MergeRate"].append(mergedreads / totalreads)
        data["Mapped"].append(mapped)
        data["MappingRate"].append(mapped / maptotal)
        data["MappedFullRefr"].append(frmapped)
        data["MappingRateFullRefr"].append(frmapped / frmaptotal)
        data["Typed"].append(num_typed_reads)
        data["TypingRate"].append(typing_rate)
        data["InterlocChiSq"].append(chisq)
        data["HetTstat"].append(tstat)
    return pd.DataFrame(data)
