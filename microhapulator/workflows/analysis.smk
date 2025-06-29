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

from importlib.resources import files
from microhapulator import api as mhapi
from microhapulator.pipe import OverviewReporter, DetailReporter
from microhapulator.profile import TypingResult
from microhapulator.thresholds import ThresholdIndex
import shutil


include: "preproc-paired.smk" if config["paired"] else "preproc-single.smk"


rule report:
    input:
        "report/img/read-mapping-qc.png",
        preproc_aux_files,
        expand("analysis/{sample}/03typing/{sample}-type.json", sample=config["samples"]),
        expand(
            "analysis/{sample}/04profiles/{sample}-{suffix}.csv",
            sample=config["samples"],
            suffix=("qual", "quant", "qual-ref", "quant-ref"),
        ),
        expand(
            "analysis/{sample}/01preprocessing/{sample}-ambig-read-counts.txt",
            sample=config["samples"],
        ),
        expand(
            "analysis/{sample}/01preprocessing/{sample}-length-filtered-read-counts.txt",
            sample=config["samples"],
        ),
        expand(
            "analysis/{sample}/02alignment/{sample}-repetitive-reads.csv",
            sample=config["samples"],
        ),
        expand(
            "analysis/{sample}/02alignment/{sample}-interlocus-balance.png",
            sample=config["samples"],
        ),
        expand(
            "analysis/{sample}/03typing/{sample}-heterozygote-balance.png",
            sample=config["samples"],
        ),
        expand("analysis/{sample}/03typing/callplots/.done", sample=config["samples"]),
        expand("analysis/{sample}/03typing/{sample}-typing-rate.tsv", sample=config["samples"]),
        expand("analysis/{sample}/03typing/{sample}-discard-rate.tsv", sample=config["samples"]),
        expand("analysis/{sample}/03typing/{sample}-gapped-rate.tsv", sample=config["samples"]),
        files("microhapulator") / "data" / "template.html",
        files("microhapulator") / "data" / "marker_details_template.html",
        files("microhapulator") / "data" / "third-party" / "bootstrap.min.css",
        files("microhapulator") / "data" / "third-party" / "fancyTable.js",
        files("microhapulator") / "data" / "third-party" / "jquery-ui.min.js",
        files("microhapulator") / "data" / "third-party" / "jquery.min.js",
    output:
        main="report/report.html",
        detail="report/marker-detail-report.html",
    run:
        thresholds = ThresholdIndex.load(
            configfile=config["thresh_file"],
            global_static=config["thresh_static"],
            global_dynamic=config["thresh_dynamic"],
            ambiguous=config["ambiguous_thresh"],
            min_read_length=config["length_thresh"],
            discard_alert=config["thresh_discard_alert"],
            gap_alert=config["thresh_gap_alert"],
        )
        reporter = OverviewReporter(
            config["samples"], thresholds, reads_are_paired=config["paired"]
        )
        with open(output.main, "w") as fh:
            print(reporter.render(), file=fh, end="")
        reporter = DetailReporter(config["samples"])
        with open(output.detail, "w") as fh:
            print(reporter.render(), file=fh, end="")
        tpdir = files("microhapulator") / "data" / "third-party"
        shell("mkdir -p report/img/")
        shell("cp -r {tpdir} report/assets")
        for sample in config["samples"]:
            shell("cp analysis/{sample}/*/*.png report/img/")
            shell("cp -r analysis/{sample}/03typing/callplots/ report/img/{sample}-callplots/")


rule copy_and_index_marker_data:
    input:
        tsv=config["mhdefn"],
        fasta=config["mhrefr"],
    output:
        fasta="marker-refr.fasta",
        index="marker-refr.mmi",
        tsv="marker-definitions.tsv",
    shell:
        """
        cp {input.tsv} {output.tsv}
        cp {input.fasta} {output.fasta}
        minimap2 -d {output.index} -k 21 -w 11 {output.fasta}
        """


rule map_sort_and_index:
    input:
        fastq="analysis/{sample}/01preprocessing/{sample}-preprocessed-reads.fastq.gz",
        index="marker-refr.mmi",
    output:
        bam="analysis/{sample}/02alignment/{sample}.bam",
        bai="analysis/{sample}/02alignment/{sample}.bam.bai",
        stats="analysis/{sample}/02alignment/{sample}.bam.stats",
    threads: 32
    shell:
        """
        minimap2 -ax sr -t {threads} {input.index} {input.fastq} | samtools view -b | samtools sort -o {output.bam}
        samtools index {output.bam}
        samtools stats {output.bam} > {output.stats}
        """


rule map_full_reference:
    input:
        index=config["hg38index"],
        fastq="analysis/{sample}/01preprocessing/{sample}-preprocessed-reads.fastq.gz",
    output:
        bam="analysis/{sample}/02alignment/{sample}-fullrefr.bam",
        bai="analysis/{sample}/02alignment/{sample}-fullrefr.bam.bai",
    threads: 32
    shell:
        """
        minimap2 -ax sr -t {threads} {input.index} {input.fastq} | samtools view -b | samtools sort -o {output.bam}
        samtools index {output.bam}
        """


rule count_fullrefr_mapped_reads:
    input:
        bam=rules.map_full_reference.output.bam,
        bai=rules.map_full_reference.output.bai,
    output:
        counts="analysis/{sample}/02alignment/{sample}-fullrefr-mapped-reads.txt",
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
        typing_result="analysis/{sample}/03typing/{sample}-type-raw.json",
    shell:
        "mhpl8r type {input} --out {output}"


rule check_typing_rate:
    input:
        json=rules.call_haplotypes.output.typing_result,
    output:
        tsv="analysis/{sample}/03typing/{sample}-typing-rate.tsv",
    run:
        result = TypingResult(fromfile=input.json)
        rates = result.typing_rate()
        rates.to_csv(output.tsv, sep="\t", index=False, float_format="%.4f")


rule check_discard_rate:
    input:
        json=rules.call_haplotypes.output.typing_result,
    output:
        tsv="analysis/{sample}/03typing/{sample}-discard-rate.tsv",
    params:
        threshold=float(config["thresh_discard_alert"]),
    run:
        result = TypingResult(fromfile=input.json)
        rates = result.typing_rate()
        rates["DiscardRate"] = 1.0 - rates.TypingRate
        rates = rates[rates.DiscardRate > params.threshold]
        rates.to_csv(output.tsv, sep="\t", index=False, float_format="%.4f")


rule check_gapped_rate:
    input:
        json=rules.call_haplotypes.output.typing_result,
    output:
        tsv="analysis/{sample}/03typing/{sample}-gapped-rate.tsv",
    params:
        threshold=float(config["thresh_gap_alert"]),
    run:
        result = TypingResult(fromfile=input.json)
        rates = result.gap_rate()
        rates = rates[rates.GappedRate > params.threshold]
        rates.to_csv(output.tsv, sep="\t", index=False, float_format="%.4f")


rule apply_filters:
    input:
        json=rules.call_haplotypes.output.typing_result,
    output:
        genotype_call="analysis/{sample}/03typing/{sample}-type.json",
    params:
        static=f"--static {config['thresh_static']}",
        dynamic=f"--dynamic {config['thresh_dynamic']}",
        threshfile=(
            "" if config["thresh_file"] in ("", None) else f"--config {config['thresh_file']}"
        ),
    shell:
        "mhpl8r filter {input} --out {output} {params.static} {params.dynamic} {params.threshfile}"


rule profile_convert:
    input:
        json=rules.apply_filters.output.genotype_call,
    output:
        qual="analysis/{sample}/04profiles/{sample}-qual.csv",
        quant="analysis/{sample}/04profiles/{sample}-quant.csv",
        qualref="analysis/{sample}/04profiles/{sample}-qual-ref.csv",
        quantref="analysis/{sample}/04profiles/{sample}-quant-ref.csv",
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
        counts="analysis/{sample}/02alignment/{sample}-marker-read-counts.csv",
        plot="analysis/{sample}/02alignment/{sample}-interlocus-balance.png",
        log="analysis/{sample}/02alignment/{sample}-interlocus-balance-chisq.txt",
    shell:
        """
        mhpl8r locbalance {input} --csv {output.counts} --figure {output.plot} --title {wildcards.sample} --quiet \
            | tee {output.log}
        """


rule heterozygote_balance:
    input:
        result=rules.apply_filters.output.genotype_call,
    output:
        counts="analysis/{sample}/03typing/{sample}-het-marker-allele-counts.csv",
        plot="analysis/{sample}/03typing/{sample}-heterozygote-balance.png",
        log="analysis/{sample}/03typing/{sample}-heterozygote-balance-pttest.txt",
    shell:
        """
        mhpl8r hetbalance {input} --csv {output.counts} --figure {output.plot} --title {wildcards.sample} --labels \
            | tee {output.log}
        """


rule plot_haplotype_calls:
    input:
        result=rules.apply_filters.output.genotype_call,
    output:
        sentinel=touch("analysis/{sample}/03typing/callplots/.done"),
    run:
        result = TypingResult(fromfile=input.result)
        mhapi.plot_haplotype_calls(
            result,
            f"analysis/{wildcards.sample}/03typing/callplots",
            sample=wildcards.sample,
        )


rule download_and_index_full_reference:
    output:
        fasta=config["hg38path"],
        index=config["hg38index"],
    shell:
        """
        mhpl8r getrefr
        """


rule repetitive_mapping:
    input:
        marker_bam="analysis/{sample}/02alignment/{sample}.bam",
        fullref_bam="analysis/{sample}/02alignment/{sample}-fullrefr.bam",
        marker_def="marker-definitions.tsv",
    output:
        counts="analysis/{sample}/02alignment/{sample}-repetitive-reads.csv",
    shell:
        """
        mhpl8r repetitive {input} --out {output}
        """


rule read_mapping_qc:
    input:
        marker=rules.map_sort_and_index.output.stats,
        full_refr=rules.count_fullrefr_mapped_reads.output.counts,
        repetitive=rules.repetitive_mapping.output.counts,
    output:
        counts="analysis/{sample}/02alignment/{sample}-read-mapping-qc.csv",
    run:
        results = mhapi.read_mapping_qc(input.marker, input.full_refr, input.repetitive)
        results.to_csv(output.counts, index=False)


rule aggregate_read_mapping_qc:
    input:
        expand(
            "analysis/{sample}/02alignment/{sample}-read-mapping-qc.csv", sample=config["samples"]
        ),
    output:
        plot="report/img/read-mapping-qc.png",
    run:
        mhapi.aggregate_read_mapping_qc(config["samples"], output.plot)
