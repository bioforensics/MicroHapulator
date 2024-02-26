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
from microhapulator.marker import MicrohapIndex
from microhapulator.pipeaux import (
    full_reference_index_files,
    final_html_report,
    aggregate_summary,
    marker_detail_report,
)
from microhapulator.profile import TypingResult
import pandas as pd
from pkg_resources import resource_filename
import shutil


include: "preproc-paired.smk" if config["paired"] else "preproc-single.smk"


rule report:
    input:
        "analysis/summary.tsv",
        preproc_aux_files,
        expand("analysis/{sample}/{sample}-type.json", sample=config["samples"]),
        expand(
            "analysis/{sample}/profiles/{sample}-{suffix}.csv",
            sample=config["samples"],
            suffix=("qual", "quant", "qual-ref", "quant-ref"),
        ),
        expand("analysis/{sample}/{sample}-interlocus-balance.png", sample=config["samples"]),
        expand("analysis/{sample}/{sample}-heterozygote-balance.png", sample=config["samples"]),
        expand("analysis/{sample}/callplots/.done", sample=config["samples"]),
        expand("analysis/{sample}/{sample}-repetitive-reads.csv", sample=config["samples"]),
        expand("analysis/{sample}/{sample}-donut.png", sample=config["samples"]),
        resource_filename("microhapulator", "data/template.html"),
        resource_filename("microhapulator", "data/marker_details_template.html"),
        resource_filename("microhapulator", "data/fancyTable.js"),
    output:
        "report.html",
        "marker-detail-report.html",
    run:
        summary = pd.read_csv("analysis/summary.tsv", sep="\t")
        final_html_report(
            config["samples"],
            summary,
            reads_are_paired=config["paired"],
            thresh_static=config["thresh_static"],
            thresh_dynamic=config["thresh_dynamic"],
            thresh_file=config["thresh_file"],
            ambiguous_read_threshold=config["ambiguous_thresh"],
        )
        marker_detail_report(config["samples"])
        jsfile = resource_filename("microhapulator", "data/fancyTable.js")
        shutil.copy(jsfile, "fancyTable.js")


rule summary:
    input:
        summary_aux_files,
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
        summary = aggregate_summary(config["samples"], config["paired"])
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
    threads: 32
    shell:
        """
        bwa mem -t {threads} {params.refr} {input.fastq} | samtools view -b | samtools sort -o {output.bam}
        samtools index {output.bam}
        """


rule count_fullrefr_mapped_reads:
    input:
        bam=rules.map_full_reference.output.bam,
        bai=rules.map_full_reference.output.bai,
    output:
        counts="analysis/{sample}/fullrefr/{sample}-fullrefr-mapped-reads.txt",
    shell:
        """
        echo -n "Total reads: " > {output.counts}
        samtools view -c -F 2304 {input.bam} >> {output.counts}
        echo -n "Mapped reads: " >> {output.counts}
        samtools view -c -F 2308 {input.bam} >> {output.counts}
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
        static=f"--static {config['thresh_static']}",
        dynamic=f"--dynamic {config['thresh_dynamic']}",
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


rule repetitive_mapping:
    input:
        marker_bam="analysis/{sample}/{sample}.bam",
        fullref_bam="analysis/{sample}/fullrefr/{sample}-fullrefr.bam",
        marker_def="marker-definitions.tsv",
    output:
        counts="analysis/{sample}/{sample}-repetitive-reads.csv",
    run:
        microhaps = MicrohapIndex.from_files(input.marker_def)
        if microhaps.has_chrom_offsets:
            shell("mhpl8r repetitive {input} --out {output}")
        else:
            shell("touch {output}")


rule read_mapping_qc:
    input:
        marker=rules.map_sort_and_index.output.counts,
        full_refr=rules.count_fullrefr_mapped_reads.output.counts,
        repetitive=rules.repetitive_mapping.output.counts,
    output:
        counts="analysis/{sample}/{sample}-read-mapping-qc.csv",
        plot="analysis/{sample}/{sample}-donut.png",
    shell:
        """
        mhpl8r mappingqc  --marker {input.marker} --refr {input.full_refr} --rep {input.repetitive}  --csv {output.counts}  --figure {output.plot} --title {wildcards.sample}
        """
