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
from microhapulator.filter import AmbigSingleReadFilter
from microhapulator.pipeaux import full_reference_index_files
from os import symlink

preproc_aux_files = chain(
    expand("analysis/{sample}/fastqc/report.html", sample=config["samples"]),
    ["analysis/read-lengths.png"],
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


rule filter_ambiguous:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        filtered="analysis/{sample}/{sample}-ambig-filtered.fastq",
        copied_fq="analysis/{sample}/preprocessed-reads.fastq",
        counts="analysis/{sample}/{sample}-ambig-read-counts.txt",
    params:
        ambig_thresh=config["ambiguous_thresh"],
    run:
        assert len(input) == 1
        out_prefix = f"analysis/{wildcards.sample}/{wildcards.sample}"
        ambig_filter = AmbigSingleReadFilter(input[0], out_prefix, params.ambig_thresh)
        ambig_filter.filter()
        ambig_filter.write_counts_output()
        shell("cp {output.filtered} {output.copied_fq}")


rule calculate_read_lengths:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        json="analysis/{sample}/{sample}-read-lengths.json",
    run:
        mhapi.calculate_read_lengths(
            input[0],
            output.json,
        )


rule plot_read_length_distributions:
    input:
        reads=expand("analysis/{sample}/{sample}-read-lengths.json", sample=config["samples"]),
    output:
        png="analysis/read-lengths.png",
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
