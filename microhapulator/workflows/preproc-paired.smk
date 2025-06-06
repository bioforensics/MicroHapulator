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

from microhapulator import api as mhapi
from microhapulator.pipe.filter import AmbigPairedReadFilter, LengthSingleReadFilter
from os import symlink
from pathlib import Path

preproc_aux_files = [
    "report/img/r1-read-lengths.png",
    "report/img/r2-read-lengths.png",
    "report/img/merged-read-lengths.png",
    "report/multiqc_report.html",
]


summary_aux_files = (expand("analysis/{sample}/flash.log", sample=config["samples"]),)


rule fastqc:
    input:
        r1="seq/{sample}_R1.fastq.gz",
        r2="seq/{sample}_R2.fastq.gz",
    output:
        r1="analysis/{sample}/01preprocessing/fastqc/R1-fastqc.html",
        r2="analysis/{sample}/01preprocessing/fastqc/R2-fastqc.html",
    params:
        outdir="analysis/{sample}/01preprocessing/fastqc/",
    threads: 2
    run:
        shell("fastqc --outdir {params.outdir} --threads {threads} {input}")
        outfiles = sorted(Path(params.outdir).glob("*.html"))
        for end, outfile in enumerate(outfiles, 1):
            outfile = Path(outfile)
            linkfile = f"{params.outdir}/R{end}-fastqc.html"
            symlink(outfile.name, linkfile)


rule multiqc:
    input:
        expand(
            "analysis/{sample}/01preprocessing/fastqc/R{end}-fastqc.html",
            sample=config["samples"],
            end=(1, 2),
        ),
    output:
        report="report/multiqc_report.html",
    shell:
        """
        multiqc analysis/*/01preprocessing/fastqc --outdir report
        """


rule filter_ambiguous:
    input:
        r1="seq/{sample}_R1.fastq.gz",
        r2="seq/{sample}_R2.fastq.gz",
    output:
        filtered_r1="analysis/{sample}/01preprocessing/{sample}-ambig-filtered-R1.fastq.gz",
        filtered_r2="analysis/{sample}/01preprocessing/{sample}-ambig-filtered-R2.fastq.gz",
        mates_r1="analysis/{sample}/01preprocessing/{sample}-ambig-R1-mates.fastq.gz",
        mates_r2="analysis/{sample}/01preprocessing/{sample}-ambig-R2-mates.fastq.gz",
        counts="analysis/{sample}/01preprocessing/{sample}-ambig-read-counts.txt",
    params:
        ambig_thresh=config["ambiguous_thresh"],
        out_prefix="analysis/{sample}/01preprocessing/{sample}",
    run:
        ambig_filter = AmbigPairedReadFilter(*input, params.out_prefix, params.ambig_thresh)
        ambig_filter.filter()
        with open(output.counts, "w") as fh:
            print(ambig_filter.summary, file=fh)


rule merge:
    input:
        r1=rules.filter_ambiguous.output.filtered_r1,
        r2=rules.filter_ambiguous.output.filtered_r2,
    output:
        mergedfq="analysis/{sample}/01preprocessing/flash/{sample}.extendedFrags.fastq.gz",
        orphan1fq="analysis/{sample}/01preprocessing/flash/{sample}.notCombined_1.fastq.gz",
        orphan2fq="analysis/{sample}/01preprocessing/flash/{sample}.notCombined_2.fastq.gz",
    log:
        "analysis/{sample}/01preprocessing/flash/flash.log",
    params:
        prefix="analysis/{sample}/01preprocessing/flash/{sample}",
    threads: 8
    shell:
        """
        flash --min-overlap=25 --max-overlap=300 --allow-outies \
            --threads={threads} --max-mismatch-density=0.25 \
            --output-prefix {params.prefix} \
            {input} \
            2>&1 | tee {log}
        bgzip {params.prefix}.extendedFrags.fastq
        bgzip {params.prefix}.notCombined_1.fastq
        bgzip {params.prefix}.notCombined_2.fastq
        """


rule filter_length:
    input:
        mergedfq=rules.merge.output.mergedfq,
    output:
        length_filtered="analysis/{sample}/01preprocessing/{sample}-length-filtered.fastq.gz",
        linkedfq="analysis/{sample}/01preprocessing/{sample}-preprocessed-reads.fastq.gz",
        counts="analysis/{sample}/01preprocessing/{sample}-length-filtered-read-counts.txt",
    params:
        length_thresh=config["length_thresh"],
    run:
        length_filter = LengthSingleReadFilter(
            input.mergedfq, output.length_filtered, params.length_thresh
        )
        length_filter.filter()
        with open(output.counts, "w") as fh:
            print(length_filter.summary, file=fh)
        symlink(Path(output.length_filtered).name, output.linkedfq)


rule calculate_read_lengths:
    input:
        r1="seq/{sample}_R1.fastq.gz",
        r2="seq/{sample}_R2.fastq.gz",
        rm=rules.merge.output.mergedfq,
    output:
        l1="analysis/{sample}/01preprocessing/{sample}-r1-read-lengths.json",
        l2="analysis/{sample}/01preprocessing/{sample}-r2-read-lengths.json",
        merged="analysis/{sample}/01preprocessing/{sample}-merged-read-lengths.json",
    run:
        mhapi.calculate_read_lengths(input.r1, output.l1)
        mhapi.calculate_read_lengths(input.r2, output.l2)
        mhapi.calculate_read_lengths(input.rm, output.merged)


rule plot_read_length_distributions:
    input:
        r1s=expand(
            "analysis/{sample}/01preprocessing/{sample}-r1-read-lengths.json",
            sample=config["samples"],
        ),
        r2s=expand(
            "analysis/{sample}/01preprocessing/{sample}-r2-read-lengths.json",
            sample=config["samples"],
        ),
        merged=expand(
            "analysis/{sample}/01preprocessing/{sample}-merged-read-lengths.json",
            sample=config["samples"],
        ),
    output:
        r1="report/img/r1-read-lengths.png",
        r2="report/img/r2-read-lengths.png",
        merged="report/img/merged-read-lengths.png",
    run:
        mhapi.read_length_dist(
            input.r1s,
            output.r1,
            config["samples"],
            config["hspace"],
            xlabel="R1 Read Length (bp)",
            color="#e41a1c",
            edgecolor="#990000",
        )
        mhapi.read_length_dist(
            input.r2s,
            output.r2,
            config["samples"],
            config["hspace"],
            xlabel="R2 Read Length (bp)",
            color="#e41a1c",
            edgecolor="#990000",
        )
        mhapi.read_length_dist(
            input.merged,
            output.merged,
            config["samples"],
            config["hspace"],
            xlabel="Merged Read Length (bp)",
            color="#377eb8",
            edgecolor="#3355aa",
        )
