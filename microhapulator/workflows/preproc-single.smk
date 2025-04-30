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
from microhapulator.pipe.filter import AmbigSingleReadFilter, LengthSingleReadFilter
from os import symlink

preproc_aux_files = ["report/img/read-lengths.png", "report/multiqc_report.html"]


summary_aux_files = list()


rule fastqc:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        html="analysis/{sample}/01preprocessing/fastqc/report.html",
    params:
        outdir="analysis/{sample}/01preprocessing/fastqc/",
    run:
        shell("fastqc --outdir {params.outdir} {input}")
        outfiles = sorted(Path(params.outdir).glob("*.html"))
        assert len(outfiles) == 1
        outfile = Path(outfiles[0])
        symlink(outfile.name, output[0])


rule multiqc:
    input:
        expand("analysis/{sample}/01preprocessing/fastqc/report.html", sample=config["samples"]),
    output:
        report="report/multiqc_report.html",
    shell:
        """
        multiqc analysis/*/01preprocessing/fastqc --outdir report
        """


rule filter_ambiguous:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        filtered_fq="analysis/{sample}/01preprocessing/{sample}-ambig-filtered.fastq.gz",
        counts="analysis/{sample}/01preprocessing/{sample}-ambig-read-counts.txt",
    params:
        ambig_thresh=config["ambiguous_thresh"],
    run:
        assert len(input) == 1
        ambig_filter = AmbigSingleReadFilter(input[0], output.filtered_fq, params.ambig_thresh)
        ambig_filter.filter()
        with open(output.counts, "w") as fh:
            print(ambig_filter.summary, file=fh)


rule filter_length:
    input:
        fq=rules.filter_ambiguous.output.filtered_fq,
    output:
        length_filtered="analysis/{sample}/01preprocessing/{sample}-length-filtered.fastq.gz",
        linkedfq="analysis/{sample}/01preprocessing/{sample}-preprocessed-reads.fastq.gz",
        counts="analysis/{sample}/01preprocessing/{sample}-length-filtered-read-counts.txt",
    params:
        length_thresh=config["length_thresh"],
    run:
        length_filter = LengthSingleReadFilter(
            input.fq, output.length_filtered, params.length_thresh
        )
        length_filter.filter()
        with open(output.counts, "w") as fh:
            print(length_filter.summary, file=fh)
        symlink(Path(output.length_filtered).name, output.linkedfq)


rule calculate_read_lengths:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        json="analysis/{sample}/01preprocessing/{sample}-read-lengths.json",
    run:
        mhapi.calculate_read_lengths(
            input[0],
            output.json,
        )


rule plot_read_length_distributions:
    input:
        reads=expand(
            "analysis/{sample}/01preprocessing/{sample}-read-lengths.json",
            sample=config["samples"],
        ),
    output:
        png="report/img/read-lengths.png",
    run:
        mhapi.read_length_dist(
            input.reads,
            output.png,
            config["samples"],
            config["hspace"],
            xlabel="Read Length (bp)",
            color="#e41a1c",
            edgecolor="#990000",
        )
