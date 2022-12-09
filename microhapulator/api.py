# -------------------------------------------------------------------------------------------------
# Copyright (c) 2021, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


from Bio import SeqIO
from collections import namedtuple, defaultdict, Counter
import json
from math import ceil
import matplotlib
from matplotlib import pyplot as plt
from microhapulator.parsers import load_marker_definitions, cross_check_marker_ids
from microhapulator.parsers import open as mhopen
from microhapulator.profile import SimulatedProfile, TypingResult
import numpy as np
import os
import pandas as pd
from pathlib import Path
import pysam
import re
from scipy.stats import chisquare, ttest_rel
from shutil import rmtree
from string import ascii_letters, digits
from subprocess import check_call, run
import sys
from tempfile import mkdtemp, TemporaryDirectory
from warnings import warn


SimulatedRead = namedtuple("SimulatedRead", ["identifier", "sequence", "quality"])


def count_and_sort(profile, include_discarded=True):
    counts = dict(
        Marker=list(),
        ReadCount=list(),
    )
    for marker, mdata in profile.data["markers"].items():
        readcount = 0
        if include_discarded:
            readcount += mdata["num_discarded_reads"]
        for haplotype, count in mdata["typing_result"].items():
            readcount += count
        counts["Marker"].append(marker)
        counts["ReadCount"].append(readcount)
    data = pd.DataFrame(counts).sort_values(["ReadCount"], ascending=False).reset_index(drop=True)
    return data


def interlocus_balance(
    result,
    include_discarded=True,
    terminal=True,
    tofile=None,
    title=None,
    figsize=(6, 4),
    dpi=200,
    color=None,
):
    """Compute interlocus balance

    Plot interlocus balance in the terminal and/or a high-resolution graphic. Also normalize read
    counts and perform a chi-square goodness-of-fit test assuming uniform read coverage across
    markers. The reported chi-square statistic measures the extent of imbalance, and can be compared
    among samples sequenced using the same panel: the minimum value of 0 represents perfectly
    uniform coverage, while the maximum value of *D* occurs when all reads map to a single marker
    (*D* represents the degrees of freedom, or the number of markers minus 1).

    :param microhapulator.profile.TypingResult result: a typing result including haplotype counts
    :param bool included_discarded: flag indicating whether to include in each marker's total read count reads that are successfully aligned but discarded because they do not span all SNPs at the marker
    :param bool terminal: flag indicating whether to print the interlocus balance histogram to standard output; enabled by default
    :param str tofile: name of image file to which the interlocus balance histogram will be written using Matplotlib; image format is inferred from file extension; by default, no image file is generated
    :param str title: add a title (such as a sample name) to the histogram plot
    :param tuple figsize: a 2-tuple of integers indicating the dimensions of the image file to be generated
    :param int dpi: resolution (in dots per inch) of the image file to be generated
    :param str color: override histogram plot color; green by default
    :return: a tuple (S, C) where S is the chi-square statistic, and C is a table of total read counts for each marker
    :rtype: tuple(float, pandas.DataFrame)
    """
    data = count_and_sort(result, include_discarded=include_discarded)
    normalized_read_counts = [c / sum(data.ReadCount) for c in data.ReadCount]
    chisq, pval = chisquare(normalized_read_counts)
    if terminal:
        with TemporaryDirectory() as tempdir:
            tfile = os.path.join(tempdir, "data.tsv")
            data.to_csv(tfile, index=False, header=False)
            run(["termgraph", tfile])
    if tofile:
        backend = matplotlib.get_backend()
        plt.switch_backend("Agg")
        plt.figure(figsize=figsize, dpi=dpi)
        x = range(len(data))
        y = [c / 1000 for c in data.ReadCount]
        if color is None:
            color = "#4daf4a"
        plt.bar(x, y, color=color)
        plt.xticks([])
        plt.xlabel("Marker", labelpad=15, fontsize=16)
        if include_discarded:
            plt.ylabel("Reads Mapped (× 1000)", labelpad=15, fontsize=16)
        else:
            plt.ylabel("Reads Mapped and Typed (× 1000)", labelpad=15, fontsize=16)
        ax = plt.gca()
        ax.yaxis.grid(True, color="#DDDDDD")
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_color("#CCCCCC")
        ax.tick_params(left=False)
        if title is not None:
            plt.title(title, pad=25, fontsize=18)
        plt.savefig(tofile, bbox_inches="tight")
        plt.switch_backend(backend)
    return chisq, data


def genotype_counts(profile):
    counts = dict(
        Marker=list(),
        Allele1Count=list(),
        Allele2Count=list(),
        Allele1Perc=list(),
        Allele2Perc=list(),
        TotalCount=list(),
    )
    for marker, mdata in profile.data["markers"].items():
        if "genotype" not in mdata:
            warn(f"Missing genotype call for marker {marker}", UserWarning)
            continue
        genotype = {d["haplotype"] for d in mdata["genotype"]}
        if len(genotype) > 2:
            warn(f"More than two alleles found for marker {marker}", UserWarning)
            continue
        if len(genotype) < 2:
            continue
        count2, count1 = sorted([mdata["typing_result"][mhallele] for mhallele in genotype])
        totalcount = count1 + count2
        counts["Marker"].append(marker)
        counts["Allele1Count"].append(count1)
        counts["Allele2Count"].append(count2)
        counts["Allele1Perc"].append(count1 / totalcount)
        counts["Allele2Perc"].append(count2 / totalcount)
        counts["TotalCount"].append(totalcount)
    data = pd.DataFrame(counts).sort_values(["TotalCount"], ascending=False).round(decimals=4)
    data = data.drop(columns=["TotalCount"]).reset_index(drop=True)
    return data


def heterozygote_balance(
    result, tofile=None, title=None, figsize=None, dpi=200, dolabels=False, absolute=False
):
    """Compute heterozygote balance

    Compute absolute and relative abundance of major and minor alleles at heterozygote loci and plot
    abundances in a high-resolution graphic. Also perform a one-sided paired t-test and report the
    t-statistic as a measure of heterozygote imbalance.

    Graphics implementation adapted from
    https://www.pythoncharts.com/matplotlib/grouped-bar-charts-matplotlib/.

    :param microhapulator.profile.TypingResult result: a filtered typing result including haplotype counts and genotype calls
    :param str tofile: name of image file to which the interlocus balance histogram will be written using Matplotlib; image format is inferred from file extension; by default, no image file is generated
    :param str title: add a title (such as a sample name) to the histogram plot
    :param tuple figsize: a 2-tuple of integers indicating the dimensions of the image file to be generated
    :param int dpi: resolution (in dots per inch) of the image file to be generated
    :param bool dolabels: flag indicating whether marker labels and read counts should be added to the figure
    :param bool absolute: plot absolute read counts rather than relative abundances
    :return: a tuple (T, C) where T is the t-statistic, and C is the table of read counts and percentages for each heterozygous marker
    :rtype: pandas.DataFrame
    """
    data = genotype_counts(result)
    tstat, pval = ttest_rel(data.Allele1Perc, data.Allele2Perc, alternative="greater")
    if tofile:
        backend = matplotlib.get_backend()
        plt.switch_backend("Agg")
        if figsize is None:
            width = len(data) / 2 if dolabels else len(data) / 4
            width = max(8, width)
            figsize = (width, 8)
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        barwidth = 0.4
        x1 = range(len(data.Marker))
        x2 = [_x + barwidth for _x in x1]
        x = [_x + barwidth / 2 for _x in x1]
        if absolute:
            y1 = data.Allele1Count / 1000
            y2 = data.Allele2Count / 1000
            ylabel = "Read Count (× 1000)"
        else:
            y1 = data.Allele1Perc
            y2 = data.Allele2Perc
            ylabel = "Percentage of Read Count"
        ax.bar(x1, y1, width=barwidth, color="#984ea3")
        ax.bar(x2, y2, width=barwidth, color="#ff7f00")
        if dolabels:
            ax.set_xticks(x)
            ax.set_xticklabels(data.Marker, rotation=90)
        else:
            ax.set_xticks([])
        ax.yaxis.grid(True, color="#DDDDDD")
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_color("#CCCCCC")
        ax.tick_params(bottom=False, left=False)
        ax.set_xlabel("Marker", labelpad=15, fontsize=16)
        ax.set_ylabel(ylabel, labelpad=15, fontsize=16)
        if title is not None:
            ax.set_title(title, pad=25, fontsize=18)
        if dolabels:
            counts = data.Allele1Count + data.Allele2Count
            for m, height, count in zip(x, y1, counts):
                ax.text(m, height + 0.01, f"{count:,}", ha="center", va="bottom", rotation=90)
        ax.legend(["Major Allele", "Minor Allele"], loc="lower left")
        plt.savefig(tofile, bbox_inches="tight")
        plt.switch_backend(backend)
    return tstat, data


def contain(prof1, prof2):
    """Perform a simple containment test

    Calculate a simple containment statistic C, representing the number of markers whose haplotypes
    in *prof1* are a subset of the haplotypes in *prof2*. Dividing by the total number of markers
    provides a rudimentary measure of containment, the proportion of *prof1* "contained" in
    *prof2*. Note that this statistic does not accommodate allelic drop-in or drop-out.

    :param microhapulator.profile.Profile prof1: a typing result or simulated genotype
    :param microhapulator.profile.Profile prof2: a typing result or simulated genotype
    :returns: a tuple (C, T), where C is the containment statistic and T is the total number of markers
    :rtype: tuple(float, int)
    """
    total = 0
    contained = 0
    for marker in prof2.markers():
        hap1 = prof1.haplotypes(marker)
        hap2 = prof2.haplotypes(marker)
        total += len(hap2)
        contained += len(hap2 & hap1)
    return contained, total


def contrib(profile):
    r"""Estimate the minimum number of DNA contributors to a suspected mixture

    Let :math:`N_{\text{al}}` represent the largest number of haplotypes observed at any marker. We
    then estimate the minimum number of DNA contributors as follows.

    :math:`C_{\text{min}} = \left\lceil\frac{N_{\text{al}}}{2}\right\rceil`

    :param microhapulator.profile.Profile profile: a typing result or a simulated genotype
    :returns: a tuple (E, N, P), where E is the estimate for the minimum number of DNA contributors, N is the number of markers supporting the estimate, and P is the percentage of markers supporting the estimate
    :rtype: tuple(int, int, float)
    """
    num_haps_per_marker = [len(profile.haplotypes(marker)) for marker in profile.markers()]
    max_num_haps = max(num_haps_per_marker)
    max_thresh = max_num_haps - 1 if max_num_haps % 2 == 0 else max_num_haps
    max_loci = sum([1 for n in num_haps_per_marker if n >= max_thresh])
    max_perc = round(max_loci / len(num_haps_per_marker), 4)
    return ceil(max_num_haps / 2), max_loci, max_perc


def diff(prof1, prof2):
    """Compare two profiles and determine the markers at which their genotypes differ

    Note: this is a generator function.

    :param microhapulator.profile.Profile prof1: typing result or simulated profile
    :param microhapulator.profile.Profile prof2: typing result or simulated profile
    :yields: for each discordant marker, a tuple (M, X, Y), where M is the marker name, X is the set of haplotypes unique to *prof1*, and Y is the set of haplotypes unique to *prof2*
    """
    allmarkers = set(prof1.markers()).union(prof2.markers())
    for marker in sorted(allmarkers):
        haps1 = prof1.haplotypes(marker)
        haps2 = prof2.haplotypes(marker)
        diff1 = haps1 - haps2
        diff2 = haps2 - haps1
        if len(diff1) > 0 or len(diff2) > 0:
            yield marker, diff1, diff2


def dist(prof1, prof2):
    """Compute a simple Hamming distance between two profiles

    :param microhapulator.profile.Profile prof1: typing result or simulated profile
    :param microhapulator.profile.Profile prof2: typing result or simulated profile
    :returns: the number of markers with a discordant genotype between the two profiles
    :rtype: int
    """
    hammdist = 0
    for marker in set(prof1.markers()).union(prof2.markers()):
        haps1 = prof1.haplotypes(marker)
        haps2 = prof2.haplotypes(marker)
        if haps1 != haps2:
            hammdist += 1
    return hammdist


def prob(frequencies, prof1, prof2=None, erate=0.01):
    r"""Compute a profile random match probability (RMP) or an RMP-based likelihood ratio (LR) test

    The LR test, when performed, assesses the relative weight of two competing propositions.

    - :math:`H_1`: the genetic profiles *prof1* and *prof2* originated from the same individuals
    - :math:`H_2`: *prof1* and *prof2* originated from two unrelated individuals in the population

    The test statistic is computed as follows.

    .. math::

        LR = \frac{P(H_1)}{P(H_2)}

    The probability :math:`P(H_1) = \epsilon^R`, where :math:`\epsilon` is the per-marker rate of
    genotyping error (*erate*) and :math:`R` is the number of markers with discordant genotypes
    between profiles. The probability :math:`P(H_2)` is the RMP of *prof1*. Note that when there is
    a perfect match between *prof1* and *prof2*, :math:`P(H_1) = 1` and the LR statistic is simply
    the reciprocal of the RMP.

    Note that the LR test as formulated assumes that *prof1* and *prof2* are close matches. The
    error rate accommodates a small amount of allelic drop-out or drop-in but is not designed to
    accommodate profiles with substantial differences.

    :param pandas.DataFrame frequencies: table of population haplotype frequencies
    :param microhapulator.profile.Profile prof1: a typing result or simulated genotype
    :param microhapulator.profile.Profile prof2: a typing result or simulated genotype; *optional*
    :param float erate: rate of genotyping error
    :returns: the RMP of `prof1` if `prof2` is undefined, or the LR test statistic if `prof2` is defined
    :rtype: float
    """
    if prof2 is None:
        return prof1.rand_match_prob(frequencies)
    else:
        return prof1.rmp_lr_test(prof2, frequencies, erate=erate)


def calc_n_reads_from_proportions(n, totalreads, prop):
    if prop is None:
        prop = [1.0 / n for _ in range(n)]
    else:
        if len(prop) != n:
            raise ValueError("mismatch between contributor number and proportions")
    normprop = [x / sum(prop) for x in prop]
    return [int(totalreads * x) for x in normprop]


def new_signature():
    return "".join([np.random.choice(list(ascii_letters + digits)) for _ in range(7)])


def sequencing(
    profile,
    markers,
    refrseqs,
    seed=None,
    threads=1,
    numreads=500000,
    readsignature=None,
    readindex=0,
):
    tempdir = mkdtemp()
    try:
        haplofile = tempdir + "/haplo.fasta"
        with open(haplofile, "w") as fh:
            for defline, sequence in profile.haploseqs(markers, refrseqs):
                print(">", defline, "\n", sequence, sep="", file=fh)
        isscmd = [
            "iss",
            "generate",
            "--n_reads",
            str(numreads),
            "--draft",
            haplofile,
            "--model",
            "MiSeq",
            "--output",
            tempdir + "/seq",
            "--quiet",
        ]
        if seed:
            isscmd.extend(["--seed", str(seed)])
        if threads:
            isscmd.extend(["--cpus", str(threads)])
        check_call(isscmd)
        f1, f2 = tempdir + "/seq_R1.fastq", tempdir + "/seq_R2.fastq"
        with open(f1, "r") as infh1, open(f2, "r") as infh2:
            if readsignature is None:
                readsignature = new_signature()
            linebuffer = list()
            for line_r1, line_r2 in zip(infh1, infh2):
                if line_r1.startswith("@mh"):
                    readindex += 1
                    prefix = f"@{readsignature}:0:0:0:0:0:{readindex} 1:N:0:0 mh"
                    line_r1 = line_r1.replace("@mh", prefix, 1)
                    prefix = f"@{readsignature}:0:0:0:0:0:{readindex} 2:N:0:0 mh"
                    line_r2 = line_r2.replace("@mh", prefix, 1)
                    line_r1 = line_r1[:-3] + "\n"
                    line_r2 = line_r2[:-3] + "\n"
                linebuffer.append((line_r1, line_r2))
                if len(linebuffer) == 4:
                    r1 = SimulatedRead(
                        identifier=linebuffer[0][0],
                        sequence=linebuffer[1][0],
                        quality=linebuffer[3][0],
                    )
                    r2 = SimulatedRead(
                        identifier=linebuffer[0][1],
                        sequence=linebuffer[1][1],
                        quality=linebuffer[3][1],
                    )
                    yield readindex, r1, r2
                    linebuffer = list()
    finally:
        rmtree(tempdir)


def seq(
    profiles,
    markers,
    refrseqs,
    seeds=None,
    threads=1,
    totalreads=500000,
    proportions=None,
    sig=None,
):
    """Simulate paired-end Illumina MiSeq sequencing of the given profile(s)

    This generator function accepts any combination of simple (single-source) or complex
    (multi-source mixture) profiles as input. Each profile is "sequenced" separately, and then all
    reads are aggregated.

    :param list profiles: list of mock profiles
    :param pandas.DataFrame markers: marker definitions, provided as a table of SNP offsets, one row per SNP; required columns: **Marker** and **Offset**, representing the distance of the SNP from the first nucleotide in the reference sequence
    :param dict refrseqs: a dictionary with marker names as the keys and marker reference sequences as the values
    :param list seeds: optional list of random seeds, one per profile
    :param int threads: number of threads to use for each sequencing task; note that optimal performance is typically achieved with a single thread
    :param int totalreads: total number of reads to generate across all profiles
    :param list proportions: optional list of relative proportions, equal to the number of profiles; by default, each profile contributes an equal number of reads
    :param str sig: an optional signature (prefix) to apply to all simulated reads
    :yields: for each read pair, a tuple (I, X, Y) where I is a serial index, X is the forward read (R1), and Y is the reverse read (R2)
    """
    n = len(profiles)
    if seeds is None:
        seeds = [np.random.randint(1, 2**32 - 1) for _ in range(n)]
    if len(seeds) != n:
        raise ValueError("number of profiles must match number of seeds")
    numreads = calc_n_reads_from_proportions(n, totalreads, proportions)
    if 0 in numreads:
        raise ValueError("specified proportions result in 0 reads for 1 or more individuals")
    readsignature = sig if sig else new_signature()
    reads_sequenced = 0
    for profile, seed, nreads in zip(profiles, seeds, numreads):
        message = "Individual seed={seed} numreads={n}".format(seed=seed, n=nreads)
        print("[MicroHapulator::seq]", message, file=sys.stderr)
        sequencer = sequencing(
            profile,
            markers,
            refrseqs,
            seed=seed,
            threads=threads,
            numreads=nreads,
            readsignature=readsignature,
            readindex=reads_sequenced,
        )
        for data in sequencer:
            yield data
        reads_sequenced = data[0]


def sim(frequencies, seed=None):
    """Simulate a diploid genotype from the specified microhaplotype frequencies

    :param pandas.DataFrame frequencies: population haplotype frequencies
    :param int seed: seed for random number generator
    :returns: a simulated genotype profile for all markers specified in the haplotype frequencies
    :rtype: microhapulator.profile.SimulatedProfile
    """
    profile = SimulatedProfile(ploidy=2)
    if seed is None:
        seed = np.random.randint(2**32 - 1)
    profile.data["metadata"] = {
        "HaploSeed": seed,
    }
    np.random.seed(seed)
    markers = sorted(frequencies.Marker.unique())
    for haploindex in range(2):
        for marker in markers:
            haplofreqs = frequencies[frequencies.Marker == marker]
            haplotypes = list(haplofreqs.Haplotype)
            freqs = list(haplofreqs.Frequency)
            freqs = [x / sum(freqs) for x in freqs]
            sampled_haplotype = np.random.choice(haplotypes, p=freqs)
            profile.add(haploindex, marker, sampled_haplotype)
    message = f"simulated microhaplotype variation at {len(markers)} markers"
    print("[MicroHapulator::sim]", message, file=sys.stderr)
    return profile


def check_index(bamfile):
    index1 = bamfile + ".bai"
    index2 = re.sub(r".bam$", ".bai", bamfile)
    for testfile in (index1, index2):
        if os.path.isfile(testfile):
            break
    else:
        pysam.index(bamfile)


def skip_read(read):
    return read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail


def tally_haplotypes(bam, offsets, minbasequal=10, max_depth=1e6):
    totaldiscarded = 0
    for locusid in sorted(offsets):
        discarded = 0
        haplotypes = defaultdict(int)
        ht = defaultdict(dict)
        varloc = set(offsets[locusid])
        cov_pos = list()
        for column in bam.pileup(locusid, min_base_quality=minbasequal, max_depth=max_depth):
            cov_pos.append(column.n)
            if column.pos not in varloc:
                continue
            for record in column.pileups:
                aligned_base = None
                if record.is_del or record.is_refskip or skip_read(record.alignment):
                    continue
                aligned_base = record.alignment.query_sequence[record.query_position]
                ht[record.alignment.query_name][column.pos] = aligned_base
        for readname, htdict in ht.items():
            htlist = [htdict[pos] for pos in sorted(htdict)]
            if len(htlist) < len(varloc):
                discarded += 1
                continue
            htstr = ",".join(htlist)
            haplotypes[htstr] += 1
        yield locusid, cov_pos, haplotypes, discarded
        totaldiscarded += discarded
    print(
        "[MicroHapulator::type] discarded",
        totaldiscarded,
        "reads with gaps or missing data at positions of interest",
        file=sys.stderr,
    )


def type(bamfile, markertsv, minbasequal=10, max_depth=1e6):
    """Perform haplotype calling

    :param str bamfile: path of a BAM file containing NGS reads aligned to marker reference sequences and sorted
    :param str markertsv: path of a TSV file containing marker metadata, specifically the offset of each SNP for every marker in the panel
    :param int minbasequal: minimum base quality (PHRED score) to be considered reliable for haplotype calling; default is 10, corresponding to Q10, i.e., 90% probability that the base call is correct
    :param float max_depth: maximum permitted read depth
    :returns: an unfiltered catalog of haplotype counts for each marker (a *typing result*)
    :rtype: microhapulator.profile.TypingResult
    """
    check_index(bamfile)
    bam = pysam.AlignmentFile(bamfile, "rb")
    markers = load_marker_definitions(markertsv)
    offsets = defaultdict(list)
    for n, row in markers.iterrows():
        offsets[row.Marker].append(row.Offset)
    cross_check_marker_ids(bam.references, offsets.keys(), "read alignments", "marker definitions")
    haplotype_caller = tally_haplotypes(bam, offsets, minbasequal=minbasequal, max_depth=max_depth)
    result = TypingResult()
    for locusid, cov_by_pos, htcounts, ndiscarded in haplotype_caller:
        result.record_coverage(locusid, cov_by_pos, ndiscarded=ndiscarded)
        for haplotype, count in htcounts.items():
            result.record_haplotype(locusid, haplotype, count)
    return result


def read_length_dist(
    fastq,
    plotfile,
    lengthsfile=None,
    xlabel="Read Length (bp)",
    xlim=None,
    scale=1000,
    title=None,
    color=None,
    edgecolor=None,
):
    """Plot distribution of read lengths

    :param str fastq: path of a FASTQ file containing NGS reads
    :param str plotfile: path of a graphic file to create
    :param str lengthsfile: if specified, read lengths will be written to the specified file
    :param str xlabel: label for the X axis
    :param tuple xlim: a 2-tuple of numbers (x1, x2) representing the start and end points of the portion of the X axis to be displayed; by default this is determined automatically
    :param float scale: scaling factor for the Y axis
    :param str title: title for the plot
    :param str color: override histogram plot color; red by default
    :param str edgecolor: override histogram edge color; dark red by default
    """
    backend = matplotlib.get_backend()
    plt.switch_backend("Agg")
    lengths = list()
    with mhopen(fastq, "r") as infh:
        for record in SeqIO.parse(infh, "fastq"):
            lengths.append(len(record))
        if lengthsfile:
            with mhopen(lengthsfile, "w") as outfh:
                json.dump(lengths, outfh)
    fig = plt.figure(figsize=(6, 4), dpi=200)
    if color is None:
        color = "#e41a1c"
    if edgecolor is None:
        edgecolor = "#990000"
    weights = [1 / scale] * len(lengths)
    plt.hist(lengths, bins=25, weights=weights, color=color, edgecolor=edgecolor)
    if xlim is None:
        xlim = (min(lengths) * 0.9, max(lengths) * 1.1)
    plt.xlim(*xlim)
    ax = plt.gca()
    ax.yaxis.grid(True, color="#DDDDDD")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#CCCCCC")
    ax.tick_params(left=False)
    ax.set_xlabel(xlabel, labelpad=15, fontsize=16)
    ax.set_ylabel(f"Frequency (× {scale})", labelpad=15, fontsize=16)
    if title:
        ax.set_title(title, pad=25, fontsize=18)
    plt.savefig(plotfile, bbox_inches="tight")
    plt.switch_backend(backend)


def plot_haplotype_calls(result, outdir, sample=None, plot_marker_name=True, ignore_low=True):
    """Plot haplotype calls for each marker in a typing result

    :param microhapulator.profile.TypingResult result: a typing result
    :param Path outdir: Path object or string indicating the directory to which the graphics files will be saved
    :param str sample: name of the sample to be included as the plot title; by default no sample name is shown
    :param boolean plot_marker_name: flag indicating whether to plot the marker name as subtitle
    :param boolean ignore_low: flag indicating whether to exclude haplotypes that fall below the detection threshold (when provided)
    """
    backend = matplotlib.get_backend()
    plt.switch_backend("Agg")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for marker in result.markers():
        count_dict = result.allele_counts(marker)
        mdata = result.data["markers"][marker]
        if ignore_low and "thresholds" in mdata and "static" in mdata["thresholds"]:
            static = mdata["thresholds"]["static"]
            count_dict = {allele: count for allele, count in count_dict.items() if count >= static}
        counts = count_dict.values()
        alleles = count_dict.keys()
        fig = plt.figure(figsize=(4, 4), dpi=150)
        plt.bar(range(len(counts)), counts, color="#999999")
        if len(counts) == 1:
            plt.xticks(range(len(counts)), labels=alleles)
        else:
            plt.xticks(range(len(counts)), labels=alleles, rotation=45, ha="right")
        plt.xlabel("Observed MH Alleles", fontsize=14)
        plt.ylabel("Read Count", fontsize=14)
        if "thresholds" in mdata and "dynamic" in mdata["thresholds"]:
            dynamic = mdata["thresholds"]["dynamic"]
            plt.axhline(y=dynamic, color="#e41a1c", linestyle="--", label=f"Dynamic Threshold")
            plt.legend()
        plt.gca().yaxis.grid(True, color="#DDDDDD")
        plt.gca().set_axisbelow(True)
        plt.gca().spines["top"].set_visible(False)
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["left"].set_visible(False)
        plt.gca().spines["bottom"].set_color("#CCCCCC")
        if plot_marker_name:
            plt.title(marker)
        if sample:
            plt.suptitle(sample)
        filename = outdir / f"{marker}.png"
        plt.savefig(filename, bbox_inches="tight")
        plt.close()
    plt.switch_backend(backend)


def get_reads_in_marker_loci(fullref_bam_file, all_marker_defs):
    fullref_bam = pysam.AlignmentFile(fullref_bam_file, "rb")
    reads_to_markers = defaultdict(list)
    for marker in set(all_marker_defs.index):
        marker_def = all_marker_defs.loc[marker]
        start = min(marker_def["OffsetHg38"])
        end = max(marker_def["OffsetHg38"]) + 1
        for read in fullref_bam.fetch(marker_def["Chrom"][0], start, end):
            if not skip_read(read):
                reads_to_markers[read.query_name].append(marker)
    return reads_to_markers


def count_repetitive_reads(marker_bam, marker, marker_def, reads_to_markers, minbasequal):
    repetitive_count = 0
    for read in marker_bam.fetch(marker):
        if not skip_read(read):
            qual_mask = np.array(read.query_qualities) >= minbasequal
            qual_filtered_positions = np.array(read.get_reference_positions(full_length=True))[
                qual_mask
            ]
            contains_varloc = bool(set(qual_filtered_positions) & set(marker_def["Offset"]))
            if contains_varloc:
                repetitive_count += marker not in reads_to_markers[read.query_name]
    return repetitive_count


def repetitive_mapping(marker_bam_file, fullref_bam_file, markertsv, minbasequal=10):
    """Count reads mapped that map to a marker sequence but preferentially align elsewhere in the full reference

    :param str marker_bam_file: path of BAM file containing read alignments to marker sequences
    :param str fullref_bam_file: path of BAM file containing read alignments to the full reference genome
    :param str markertsv: path of a TSV file containing marker metadata including the offset of each SNP for every marker in the panel and the chromosome and coordinate of each in the reference genome
    :param int minbasequal: minimum base quality (PHRED score) to be considered reliable for haplotype calling; default is 10, corresponding to Q10, i.e., 90% probability that the base call is correct
    """
    counts = dict(
        Marker=list(),
        RepetitiveReads=list(),
    )
    marker_bam = pysam.AlignmentFile(marker_bam_file, "rb")
    all_marker_defs = load_marker_definitions(markertsv).set_index("Marker")
    reads_to_marker = get_reads_in_marker_loci(fullref_bam_file, all_marker_defs)
    for marker in set(all_marker_defs.index):
        marker_def = all_marker_defs.loc[marker]
        repetitive_count = count_repetitive_reads(
            marker_bam, marker, marker_def, reads_to_marker, minbasequal=minbasequal
        )
        counts["Marker"].append(marker)
        counts["RepetitiveReads"].append(repetitive_count)
    data = pd.DataFrame(counts).sort_values(by="Marker").reset_index(drop=True)
    return data


def read_mapping_qc(marker_mapped, refr_mapped, repetitive_mapped, figure, sample=None):
    """Count on target, off target, repetitive, and contaminant reads
    :param str marker_mapped: path of txt file containing number of reads mapped to marker sequences
    :param str refr_mapped: path of txt file containing number of reads mapped to the full human reference genome
    :param str repetitive_mapped: path of txt file containing number of reads mapped preferentially to non-marker loci in the human reference genome
    :param str output: path where the png file of the plot will be saved
    :param str sample: name of the sample to be included as the plot title; by default no sample name is shown
    
    """
    data=count_mapped_read_types(marker_mapped, refr_mapped, repetitive_mapped)
    percs =  data.values[0] / sum(data.values[0]) * 100
    backend = matplotlib.get_backend()
    plt.switch_backend("Agg")
    plt.figure(figsize=(18,6))
    labels = ["on target",  "off target", "contamination", "repetitive"]
    wedges, texts = plt.pie(percs, startangle=1)
    kw = dict(arrowprops=dict(arrowstyle="-", linewidth=1), zorder=0, va="center")
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2.0 + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        plt.annotate(
            f'{labels[i]} ({"{:.2f}".format(percs[i])}%)',
            xy=(x, y),
            xytext=(1.4 * np.sign(x), 1.4 * y),
            horizontalalignment=horizontalalignment,
            fontsize=24,
            **kw,
        )
    circle = plt.Circle((0, 0), 0.7, color="white")
    plt.gca().add_artist(circle)
    plt.title(sample, fontsize=24, fontweight="bold")
    plt.tight_layout()
    plt.savefig(figure, dpi=300)
    plt.switch_backend(backend)
    return data

    
def count_mapped_read_types(marker_mapped, refr_mapped, repetitive_mapped):
    counts = list()
    with open(marker_mapped) as fh1:
        total_reads = int(next(fh1).strip().split(": ")[1])
        marker_mapped_count = int(next(fh1).strip().split(": ")[1])
    with open(refr_mapped) as fh2:
        next(fh2)
        refr_mapped_count = int(next(fh2).strip().split(": ")[1])
    try:
        repetitive_df = pd.read_csv(repetitive_mapped)
        repetitive_count = np.sum(repetitive_df["RepetitiveReads"])
        on_target_count = marker_mapped_count - repetitive_count
    except Exception:
        repetitive_count = None
        on_target_count = marker_mapped_count
    contam_count = total_reads - refr_mapped_count
    off_target_count = refr_mapped_count - marker_mapped_count
    counts.append(on_target_count)
    counts.append(off_target_count)
    counts.append(contam_count)
    if repetitive_count:
        counts.append(repetitive_count)
    cols=["OnTarget", "OffTarget", "Contaminant", "Repetitive"]
    return pd.DataFrame([counts], columns=cols[:len(counts)]) 
