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


rule filter_ambiguous:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
    output:
        filtered="analysis/{sample}/{sample}-ambig-filtered.fastq",
        copied_fq="analysis/{sample}/reads.fastq",
        counts="analysis/{sample}/{sample}-ambig-read-counts.txt",
    params:
        ambig_thresh=config["ambiguous_thresh"],
    run:
        assert len(input) == 1
        out_prefix = f"analysis/{wildcards.sample}/{wildcards.sample}"
        ambig_filter = AmbigSingleReadFilter(input[0], out_prefix, params.ambig_thresh)
        ambig_filter.filter()
        shell("cp {output.filtered} {output.copied_fq}")


rule read_length_distributions:
    input:
        lambda wildcards: sorted([fq for fq in config["readfiles"] if wildcards.sample in fq]),
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
