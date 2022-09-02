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
from microhapulator.parsers import cross_check_marker_ids
from microhapulator.parsers import load_marker_reference_sequences, load_marker_definitions
from os import cpu_count, symlink
from pathlib import Path
from pkg_resources import resource_filename
from shutil import copy
from snakemake import snakemake
import sys


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


def validate_sample_input_files(numfiles, sample, paired_end_ok=True, single_end_ok=True):
    if not paired_end_ok and not single_end_ok:
        raise ValueError(f"both paired-end and single-end data disabled, cannot proceed")
    if numfiles == 0:
        raise FileNotFoundError(f"sample {sample}: found 0 FASTQ files")
    elif numfiles not in (1, 2):
        raise ValueError(f"sample {sample}: found {numfiles} FASTQ files, not supported")
    elif numfiles == 2 and not paired_end_ok:
        raise ValueError(f"sample {sample}: found 2 FASTQ files, but paired-end mode disabled")
    elif numfiles == 1 and not single_end_ok:
        raise ValueError(f"sample {sample}: found 1 FASTQ file, but single-end mode disabled")
    return True


def get_input_files(sample_names, seqpath, paired_end_ok=True, single_end_ok=True, suffixes=None):
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
        filelist_ok = validate_sample_input_files(
            len(filelist), sample, paired_end_ok=paired_end_ok, single_end_ok=single_end_ok
        )
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
    # Normally the defaults for the next three arguments would be `None`. But the values are passed
    # to Snakemake, which intermittently casts NoneTypes as a string, causing issues with the
    # downstream code; e.g., attempts to open a non-existent `None` file. The empty string defaults
    # and the lambda functions for data types are a workaround.
    cli.add_argument(
        "-s",
        "--static",
        metavar="ST",
        type=lambda a: "" if a == "" else int(a),
        default="",
        help="global fixed read count threshold",
    )
    cli.add_argument(
        "-d",
        "--dynamic",
        metavar="DT",
        type=lambda a: "" if a == "" else float(a),
        default="",
        help="global percentage of total read count threshold; e.g. use --dynamic=0.02 to apply a 2%% analytical threshold",
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
        dest="paired_end_ok",
        action="store_false",
        help="accept single-end reads only; by default, single-end or paired-end is accepted on a per-sample basis",
    )
    cli.add_argument(
        "--paired",
        dest="single_end_ok",
        action="store_false",
        help="accept paired-end reads only; by default, single-end or paired-end is accepted on a per-sample basis",
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
    defn = load_marker_definitions(markerdefn)
    seqs = load_marker_reference_sequences(markerseqs)
    print("[MicroHapulator] validating panel configuration", file=sys.stderr)
    cross_check_marker_ids(
        seqs.keys(), defn.Marker, "marker reference sequences", "marker definitions"
    )


def main(args):
    validate_panel_config(args.markerrefr, args.markerdefn)
    samples, filenames = get_input_files(
        args.samples, args.seqpath, args.paired_end_ok, args.single_end_ok
    )
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
    )
    snakefile = resource_filename("microhapulator", "Snakefile")
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
