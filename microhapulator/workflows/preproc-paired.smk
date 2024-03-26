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

from glob import glob
from itertools import chain
from microhapulator import api as mhapi
from microhapulator.filter import AmbigPairedReadFilter
from microhapulator.pipeaux import full_reference_index_files
from os import symlink

preproc_aux_files = chain(
    expand("analysis/{sample}/fastqc/R{end}-fastqc.html", sample=config["samples"], end=(1, 2)),
    [
        "analysis/r1-read-lengths.png",
        "analysis/r2-read-lengths.png",
        "analysis/merged-read-lengths.png",
    ],
)

summary_aux_files = (expand("analysis/{sample}/flash.log", sample=config["samples"]),)


rule fastqc:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        r1="analysis/{sample}/fastqc/R1-fastqc.html",
        r2="analysis/{sample}/fastqc/R2-fastqc.html",
    threads: 2
    run:
        shell("fastqc --outdir analysis/{wildcards.sample}/fastqc/ --threads {threads} {input}")
        outfiles = sorted(glob(f"analysis/{wildcards.sample}/fastqc/*.html"))
        for end, outfile in enumerate(outfiles, 1):
            outfile = Path(outfile)
            linkfile = f"analysis/{wildcards.sample}/fastqc/R{end}-fastqc.html"
            symlink(outfile.name, linkfile)


rule filter_ambiguous:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        filtered_r1="analysis/{sample}/{sample}-ambig-filtered-R1.fastq",
        filtered_r2="analysis/{sample}/{sample}-ambig-filtered-R2.fastq",
        mates_r1="analysis/{sample}/{sample}-ambig-R1-mates.fastq",
        mates_r2="analysis/{sample}/{sample}-ambig-R2-mates.fastq",
        counts="analysis/{sample}/{sample}-ambig-read-counts.txt",
    params:
        ambig_thresh=config["ambiguous_thresh"],
        out_prefix="analysis/{sample}/{sample}",
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
        mergedfq="analysis/{sample}/{sample}.extendedFrags.fastq",
        linkedfq="analysis/{sample}/{sample}-preprocessed-reads.fastq",
        log="analysis/{sample}/flash.log",
    threads: 8
    shell:
        """
        flash --min-overlap=25 --max-overlap=300 --allow-outies \
            --threads={threads} --max-mismatch-density=0.25 \
            --output-prefix analysis/{wildcards.sample}/{wildcards.sample} \
            {input} \
            2>&1 | tee {output.log}
        ln -s {wildcards.sample}.extendedFrags.fastq {output.linkedfq}
        """


rule calculate_read_lengths:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
        rules.merge.output.linkedfq,
    output:
        l1="analysis/{sample}/{sample}-r1-read-lengths.json",
        l2="analysis/{sample}/{sample}-r2-read-lengths.json",
        merged="analysis/{sample}/{sample}-merged-read-lengths.json",
    run:
        mhapi.calculate_read_lengths(
            input[0],
            output.l1,
        )
        mhapi.calculate_read_lengths(
            input[1],
            output.l2,
        )
        mhapi.calculate_read_lengths(input[2], output.merged)


rule plot_read_length_distributions:
    input:
        r1s=expand("analysis/{sample}/{sample}-r1-read-lengths.json", sample=config["samples"]),
        r2s=expand("analysis/{sample}/{sample}-r2-read-lengths.json", sample=config["samples"]),
        merged=expand(
            "analysis/{sample}/{sample}-merged-read-lengths.json", sample=config["samples"]
        ),
    output:
        r1="analysis/r1-read-lengths.png",
        r2="analysis/r2-read-lengths.png",
        merged="analysis/merged-read-lengths.png",
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
