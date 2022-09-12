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
from microhapulator import api as mhapi
from microhapulator.pipeaux import full_reference_index_files
from os import symlink


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


rule merge:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        mergedfq="analysis/{sample}/{sample}.extendedFrags.fastq",
        log="analysis/{sample}/flash.log",
    threads: 8
    shell:
        """
        flash --min-overlap=25 --max-overlap=300 --allow-outies \
            --threads={threads} --max-mismatch-density=0.25 \
            --output-prefix analysis/{wildcards.sample}/{wildcards.sample} \
            {input} \
            2>&1 | tee {output.log}
        """


rule read_length_distributions:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
        rules.merge.output.mergedfq,
    output:
        r1="analysis/{sample}/{sample}-r1-read-lengths.png",
        r2="analysis/{sample}/{sample}-r2-read-lengths.png",
        l1="analysis/{sample}/{sample}-r1-read-lengths.json",
        l2="analysis/{sample}/{sample}-r2-read-lengths.json",
        merged="analysis/{sample}/{sample}-merged-read-lengths.png",
    run:
        mhapi.read_length_dist(
            input[0],
            output.r1,
            lengthsfile=output.l1,
            title=f"{wildcards.sample} R1",
            color="#e41a1c",
            edgecolor="#990000",
        )
        mhapi.read_length_dist(
            input[1],
            output.r2,
            lengthsfile=output.l2,
            title=f"{wildcards.sample} R2",
            color="#e41a1c",
            edgecolor="#990000",
        )
        mhapi.read_length_dist(
            input[2],
            output.merged,
            xlabel="Length of Merged Read Pair (bp)",
            title=wildcards.sample,
            color="#377eb8",
            edgecolor="#000099",
        )


rule map_sort_and_index:
    input:
        fastq=rules.merge.output.mergedfq,
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
        fastq=rules.merge.output.mergedfq,
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
