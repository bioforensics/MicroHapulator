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

preproc_aux_files = chain(
    expand("analysis/{sample}/{sample}-read-lengths.png", sample=config["samples"]),
    expand("analysis/{sample}/fastqc/report.html", sample=config["samples"]),
)

summary_aux_files = list()


rule fastqc:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        html="analysis/{sample}/fastqc/report.html",
    run:
        shell("fastqc --outdir analysis/{wildcards.sample}/fastqc/ --threads {threads} {input}")
        outfiles = sorted(glob(f"analysis/{wildcards.sample}/fastqc/*.html"))
        assert len(outfiles) == 1
        outfile = Path(outfiles[0])
        linkfile = f"analysis/{wildcards.sample}/fastqc/report.html"
        symlink(outfile.name, linkfile)


rule fastq_reads:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        fastq="analysis/{sample}/reads.fastq",
    run:
        assert len(input) == 1
        if input[0].endswith(".gz"):
            shell("gunzip -c {input[0]} > {output.fastq}")
        else:
            shell("cp {input[0]} {output.fastq}")


rule read_length_distributions:
    input:
        rules.fastq_reads.output.fastq,
    output:
        png="analysis/{sample}/{sample}-read-lengths.png",
        json="analysis/{sample}/{sample}-read-lengths.json",
    run:
        mhapi.read_length_dist(
            input[0],
            output.png,
            lengthsfile=output.json,
            title=wildcards.sample,
            color="#e41a1c",
            edgecolor="#990000",
        )
