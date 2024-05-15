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


from argparse import SUPPRESS
from collections import defaultdict
from microhapulator.marker import MicrohapIndex
from os import cpu_count, symlink
from pathlib import Path
from pkg_resources import resource_filename
from shutil import copy
from snakemake import snakemake
import sys
from warnings import warn


def check_sample_names(samples):
    """Parse sample names and perform sanity checks

    Input is expected to be a list of strings corresponding to sample names. This function ensures
    that there are no duplicate sample names and no sample names that are substrings of other
    sample names.

    If the list contains only a single string and that string corresponds to a valid file path, it
    is assumed to be a file containing a list of sample names, one per line.
    """
    if len(samples) == 1 and Path(samples[0]).is_file():
        with open(samples[0], "r") as fh:
            samples = fh.read().strip().split("\n")
    samples = list(sorted(samples))
    for s1 in samples:
        for s2 in samples:
            if s1 == s2:
                continue
            if s1 in s2:
                message = f"cannot correctly process a sample name that is a substring of another sample name: {s1} vs. {s2}"
                raise ValueError(message)
    return samples


def traverse(dirpath):
    dirpath = Path(dirpath)
    if not dirpath.is_dir():
        return
    for subpath in dirpath.iterdir():
        if subpath.is_dir():
            yield from traverse(subpath)
        else:
            yield subpath


def validate_sample_input_files(numfiles, sample, reads_are_paired=True):
    if numfiles == 0:
        raise FileNotFoundError(f"sample {sample}: found 0 FASTQ files")
    if reads_are_paired:
        exp_num_fastq_files = 2
        mode = "paired"
    else:
        exp_num_fastq_files = 1
        mode = "single"
    if numfiles != exp_num_fastq_files:
        message = (
            f"sample {sample}: found {numfiles} FASTQ files"
            f", expected {exp_num_fastq_files} in {mode}-end mode"
        )
        raise ValueError(message)
    return True


def get_input_files(sample_names, seqpath, reads_are_paired=True, suffixes=None):
    """Find input files for each sample

    This function traverses `seqpath` and any of its sub-directories for FASTQ files. Any FASTQ
    file matching one of the sample names is stored in a dictionary with that sample name. Then the
    function checks each sample to test whether the number of files found matches the number of
    expected files. Files for samples that pass this test are stored in a list to be passed to
    Snakemake. (I would prefer to pass the data as a dictionary, but Snakemake complains when the
    config object is a nested dictionary. So instead we'll reconstruct the dictionary from this
    list in the Snakefile.)
    """
    if suffixes is None:
        suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
    sample_names = check_sample_names(sample_names)
    files = defaultdict(list)
    for filepath in traverse(seqpath):
        if not filepath.name.endswith(suffixes):
            continue
        for sample in sample_names:
            if sample in filepath.name:
                files[sample].append(filepath)
    final_file_list = list()
    for sample in sample_names:
        filelist = files[sample]
        filelist_ok = validate_sample_input_files(len(filelist), sample, reads_are_paired)
        if filelist_ok:
            final_file_list.extend(filelist)
    unique_file_names = set([filepath.name for filepath in final_file_list])
    if len(unique_file_names) != len(final_file_list):
        raise ValueError("duplicate FASTQ file names found; refusing to proceed")
    return sample_names, final_file_list


def link_or_copy_input(filenames, workdir, docopy=False):
    seqdir = Path(workdir) / "seq"
    seqdir.mkdir(parents=True, exist_ok=True)
    workingfiles = list()
    for filepath in filenames:
        workingpath = seqdir / filepath.name
        workingfiles.append(Path("seq") / workingpath.name)
        if workingpath.is_file():
            continue
        createfunc = copy if docopy else symlink
        createfunc(filepath.resolve(), workingpath)
    return [str(wf) for wf in workingfiles]


def subparser(subparsers):
    desc = "Perform a complete end-to-end microhap analysis pipeline"
    cli = subparsers.add_parser("pipe", description=desc)
    cli.add_argument(
        "-w",
        "--workdir",
        metavar="D",
        default=".",
        help="pipeline working directory; default is current directory",
    )
    cli.add_argument(
        "-n",
        "--dryrun",
        action="store_true",
        help="do not execute the workflow, but display what would have been done",
    )
    cli.add_argument(
        "-t",
        "--threads",
        metavar="T",
        type=int,
        default=cpu_count(),
        help="process each batch using T threads; by default, one thread per available core is used",
    )
    cli.add_argument(
        "-s",
        "--static",
        metavar="ST",
        type=int,
        default=5,
        help="global fixed read count threshold; ST=5 by default",
    )
    cli.add_argument(
        "-d",
        "--dynamic",
        metavar="DT",
        type=float,
        default=0.02,
        help="global percentage of total read count threshold; e.g. use --dynamic=0.02 to apply a 2%% analytical threshold; DT=0.02 by default",
    )
    cli.add_argument(
        "-a",
        "--ambiguous-thresh",
        metavar="AT",
        type=float,
        default=0.2,
        help="filter out reads with more than AT percent of ambiguous characters ('N'); AT=0.2 by default",
    )
    cli.add_argument(
        "-l",
        "--length-thresh",
        metavar="LT",
        type=float,
        default=50,
        help="filter out reads that are less than LT bp long; LT=50 by default",
    )
    cli.add_argument(
        "-c",
        "--config",
        metavar="CSV",
        default="",
        help="CSV file specifying marker-specific thresholds to override global thresholds; three required columns: 'Marker' for the marker name; 'Static' and 'Dynamic' for marker-specific thresholds",
    )
    cli.add_argument(
        "--single",
        dest="reads_are_paired",
        action="store_false",
        help="accept single-end reads only; by default, only paired-end reads are accepted",
    )
    cli.add_argument(
        "--copy-input",
        action="store_true",
        help="copy input files to working directory; by default, input files are symlinked",
    )
    cli.add_argument(
        "--hg38",
        default=resource_filename("microhapulator", "data/hg38.fasta.gz"),
        help=SUPPRESS,
        # Hidden option for testing purposes
    )
    cli.add_argument(
        "--hspace",
        metavar="HS",
        default=-0.7,
        help="horizontal spacing between samples in the read distribution length ridge plots; negative value for this parameter enables overlapping plots; HS=-0.7 by default",
    )
    cli.add_argument(
        "markerrefr", help="path to a FASTA file containing marker reference sequences"
    )
    cli.add_argument("markerdefn", help="path to a TSV file containing marker definitions")
    cli.add_argument("seqpath", help="path to a directory containing FASTQ files")
    cli.add_argument(
        "samples",
        nargs="+",
        help="list of sample names or path to .txt file containing sample names",
    )


def validate_panel_config(markerseqs, markerdefn):
    print("[MicroHapulator] validating panel configuration", file=sys.stderr)
    index = MicrohapIndex.from_files(markerdefn, markerseqs)
    index.validate()


def main(args):
    validate_panel_config(args.markerrefr, args.markerdefn)
    samples, filenames = get_input_files(args.samples, args.seqpath, args.reads_are_paired)
    workingfiles = link_or_copy_input(filenames, args.workdir, docopy=args.copy_input)
    config = dict(
        samples=samples,
        readfiles=workingfiles,
        mhrefr=Path(args.markerrefr).resolve(),
        mhdefn=Path(args.markerdefn).resolve(),
        hg38path=Path(args.hg38).resolve(),
        thresh_static=args.static,
        thresh_dynamic=args.dynamic,
        thresh_file=args.config,
        paired=args.reads_are_paired,
        ambiguous_thresh=args.ambiguous_thresh,
        length_thresh=args.length_thresh,
        hspace=args.hspace,
    )
    snakefile = resource_filename("microhapulator", "workflows/analysis.smk")
    success = snakemake(
        snakefile,
        cores=args.threads,
        printshellcmds=True,
        dryrun=args.dryrun,
        config=config,
        workdir=args.workdir,
    )
    if not success:
        return 1
