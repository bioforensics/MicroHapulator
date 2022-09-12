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

from matplotlib import pyplot as plt
from microhapulator import api as mhapi
from microhapulator.pipeaux import full_reference_index_files, final_html_report, aggregate_summary
from microhapulator.profile import TypingResult
import pandas as pd
from pkg_resources import resource_filename


include: "preproc-paired.smk" if config["paired"] else "preproc-single.smk"


rule report:
    input:
        "analysis/summary.tsv",
        preproc_files,
        expand("analysis/{sample}/{sample}-type.json", sample=config["samples"]),
        expand(
            "analysis/{sample}/profiles/{sample}-{suffix}.csv",
            sample=config["samples"],
            suffix=("qual", "quant", "qual-ref", "quant-ref"),
        ),
        expand("analysis/{sample}/{sample}-interlocus-balance.png", sample=config["samples"]),
        expand("analysis/{sample}/{sample}-heterozygote-balance.png", sample=config["samples"]),
        expand("analysis/{sample}/callplots/.done", sample=config["samples"]),
        expand("analysis/{sample}/{sample}-off-target-reads.csv", sample=config["samples"]),
        resource_filename("microhapulator", "data/template.html"),
    output:
        "report.html",
    run:
        summary = pd.read_csv("analysis/summary.tsv", sep="\t")
        final_html_report(config["samples"], summary)


rule summary:
    input:
        expand("analysis/{sample}/flash.log", sample=config["samples"]),
        expand("analysis/{sample}/{sample}-mapped-reads.txt", sample=config["samples"]),
        expand(
            "analysis/{sample}/fullrefr/{sample}-fullrefr-mapped-reads.txt",
            sample=config["samples"],
        ),
        expand("analysis/{sample}/{sample}-typing-rate.tsv", sample=config["samples"]),
        expand("analysis/{sample}/{sample}-interlocus-balance-chisq.txt", sample=config["samples"]),
        expand(
            "analysis/{sample}/{sample}-heterozygote-balance-pttest.txt",
            sample=config["samples"],
        ),
    output:
        tsv="analysis/summary.tsv",
    run:
        summary = aggregate_summary(config["samples"])
        summary.to_csv(output.tsv, sep="\t", index=False, float_format="%.4f")


rule copy_and_index_marker_data:
    input:
        tsv=config["mhdefn"],
        fasta=config["mhrefr"],
    output:
        expand("marker-refr.fasta.{suffix}", suffix=("amb", "ann", "bwt", "pac", "sa")),
        fasta="marker-refr.fasta",
        tsv="marker-definitions.tsv",
    shell:
        """
        cp {input.tsv} {output.tsv}
        cp {input.fasta} {output.fasta}
        bwa index {output.fasta}
        """


rule map_sort_and_index:
    input:
        fastq="analysis/{sample}/reads.fastq",
        fasta="marker-refr.fasta",
        idx=expand("marker-refr.fasta.{suffix}", suffix=("amb", "ann", "bwt", "pac", "sa")),
    output:
        bam="analysis/{sample}/{sample}.bam",
        bai="analysis/{sample}/{sample}.bam.bai",
        counts="analysis/{sample}/{sample}-mapped-reads.txt",
    threads: 32
    shell:
        """
        bwa mem -t {threads} {input.fasta} {input.fastq} | samtools view -b | samtools sort -o {output.bam}
        samtools index {output.bam}
        echo -n "Total reads: " > {output.counts}
        samtools view -c -F 2304 {output.bam} >> {output.counts}
        echo -n "Mapped reads: " >> {output.counts}
        samtools view -c -F 2308 {output.bam} >> {output.counts}
        """


rule map_full_reference:
    params:
        refr=config["hg38path"],
    input:
        full_reference_index_files(config["hg38path"]),
        fastq="analysis/{sample}/reads.fastq",
    output:
        bam="analysis/{sample}/fullrefr/{sample}-fullrefr.bam",
        bai="analysis/{sample}/fullrefr/{sample}-fullrefr.bam.bai",
        counts="analysis/{sample}/fullrefr/{sample}-fullrefr-mapped-reads.txt",
    threads: 32
    shell:
        """
        bwa mem -t {threads} {params.refr} {input.fastq} | samtools view -b | samtools sort -o {output.bam}
        samtools index {output.bam}
        echo -n "Total reads: " > {output.counts}
        samtools view -c -F 2304 {output.bam} >> {output.counts}
        echo -n "Mapped reads: " >> {output.counts}
        samtools view -c -F 2308 {output.bam} >> {output.counts}
        """


rule call_haplotypes:
    input:
        marker_defn=config["mhdefn"],
        bam=rules.map_sort_and_index.output.bam,
    output:
        typing_result="analysis/{sample}/{sample}-type-raw.json",
    shell:
        "mhpl8r type {input} --out {output}"


rule check_typing_rate:
    input:
        json=rules.call_haplotypes.output.typing_result,
    output:
        tsv="analysis/{sample}/{sample}-typing-rate.tsv",
    run:
        result = TypingResult(fromfile=input.json)
        rates = result.typing_rate()
        rates.to_csv(output.tsv, sep="\t", index=False, float_format="%.4f")


rule apply_filters:
    input:
        json=rules.call_haplotypes.output.typing_result,
    output:
        genotype_call="analysis/{sample}/{sample}-type.json",
    params:
        static="" if config["thresh_static"] == "" else f"--static {config['thresh_static']}",
        dynamic="" if config["thresh_dynamic"] == "" else f"--dynamic {config['thresh_dynamic']}",
        threshfile="" if config["thresh_file"] == "" else f"--config {config['thresh_file']}",
    shell:
        "mhpl8r filter {input} --out {output} {params.static} {params.dynamic} {params.threshfile}"


rule profile_convert:
    input:
        json=rules.apply_filters.output.genotype_call,
    output:
        qual="analysis/{sample}/profiles/{sample}-qual.csv",
        quant="analysis/{sample}/profiles/{sample}-quant.csv",
        qualref="analysis/{sample}/profiles/{sample}-qual-ref.csv",
        quantref="analysis/{sample}/profiles/{sample}-quant-ref.csv",
    shell:
        """
        mhpl8r convert {input} {wildcards.sample} --out {output[0]} --no-counts
        mhpl8r convert {input} {wildcards.sample} --out {output[1]}
        mhpl8r convert {input} {wildcards.sample} --out {output[2]} --no-counts --fix-homo
        mhpl8r convert {input} {wildcards.sample} --out {output[3]} --fix-homo
        """


rule interlocus_balance:
    input:
        result=rules.apply_filters.output.genotype_call,
    output:
        counts="analysis/{sample}/{sample}-marker-read-counts.csv",
        plot="analysis/{sample}/{sample}-interlocus-balance.png",
        log="analysis/{sample}/{sample}-interlocus-balance-chisq.txt",
    shell:
        """
        mhpl8r locbalance {input} --csv {output.counts} --figure {output.plot} --title {wildcards.sample} --quiet \
            | tee {output.log}
        """


rule heterozygote_balance:
    input:
        result=rules.apply_filters.output.genotype_call,
    output:
        counts="analysis/{sample}/{sample}-het-marker-allele-counts.csv",
        plot="analysis/{sample}/{sample}-heterozygote-balance.png",
        log="analysis/{sample}/{sample}-heterozygote-balance-pttest.txt",
    shell:
        """
        mhpl8r hetbalance {input} --csv {output.counts} --figure {output.plot} --title {wildcards.sample} --labels \
            | tee {output.log}
        """


rule plot_haplotype_calls:
    input:
        result=rules.apply_filters.output.genotype_call,
    output:
        sentinel=touch("analysis/{sample}/callplots/.done"),
    run:
        result = TypingResult(fromfile=input.result)
        mhapi.plot_haplotype_calls(
            result, f"analysis/{wildcards.sample}/callplots", sample=wildcards.sample
        )


rule download_and_index_full_reference:
    output:
        full_reference_index_files(config["hg38path"]),
    shell:
        """
        echo 'WARNING: if you have not previously downloaded and indexed the GRCh38 assembly for use with MicroHapulator, this can take an hour or more!!!'
        mhpl8r getrefr
        """


rule off_target_mapping:
    input:
        marker_bam=rules.map_sort_and_index.output.bam,
        fullref_bam=rules.map_full_reference.output.bam,
        marker_def=rules.copy_and_index_marker_data.output.tsv,
    output:
        counts="analysis/{sample}/{sample}-off-target-reads.csv",
    shell:
        """
        mhpl8r offtarget {input} --out {output}
        """
