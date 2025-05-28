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
import ezfastq
from ezfastq.pair import PairMode
from importlib.resources import files
from microhapulator.marker import MicrohapIndex
from os import cpu_count
from pathlib import Path
from snakemake.api import SnakemakeApi as smk_api
from snakemake.exceptions import WorkflowError
from snakemake.settings import types as smk_types
import sys


def main(args):
    validate_panel_config(args.markerrefr, args.markerdefn)
    pair_mode = PairMode.PairedEnd if args.reads_are_paired else PairMode.SingleEnd
    ezfastq.api.copy(args.samples, args.seqpath, pair_mode=pair_mode, workdir=Path(args.workdir))
    config = dict(
        samples=args.samples,
        mhrefr=str(Path(args.markerrefr).resolve()),
        mhdefn=str(Path(args.markerdefn).resolve()),
        hg38path=str(Path(args.hg38).resolve()),
        hg38index=str(Path(args.hg38idx).resolve()),
        thresh_static=args.static,
        thresh_dynamic=args.dynamic,
        thresh_file=args.config,
        paired=args.reads_are_paired,
        ambiguous_thresh=args.ambiguous_thresh,
        length_thresh=args.length_thresh,
        thresh_discard_alert=args.discard_alert,
        thresh_gap_alert=args.gap_alert,
        hspace=args.hspace,
    )
    try:
        run_snakemake(config, args.workdir, cores=args.threads, dryrun=args.dryrun)
        return 0
    except WorkflowError:
        return 1


def run_snakemake(config, workdir, cores=1, dryrun=False):
    outcfg = smk_types.OutputSettings(printshellcmds=True)
    with smk_api(outcfg) as smk:
        workflow = smk.workflow(
            config_settings=smk_types.ConfigSettings(config=config),
            storage_settings=smk_types.StorageSettings(),
            resource_settings=smk_types.ResourceSettings(cores=cores),
            snakefile=files("microhapulator") / "workflows" / "analysis.smk",
            workdir=Path(workdir),
        )
        dag = workflow.dag(smk_types.DAGSettings(targets=["report"]))
        mode = "dryrun" if dryrun else "local"
        dag.execute_workflow(executor=mode)


def validate_panel_config(markerseqs, markerdefn):
    print("[MicroHapulator] validating panel configuration", file=sys.stderr)
    index = MicrohapIndex.from_files(markerdefn, markerseqs)
    index.validate()


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
        "-D",
        "--discard-alert",
        metavar="DA",
        type=float,
        default=0.25,
        help="issue an alert in the final report for each marker whose read discard rate (proportion of reads that could not be typed) exceeds DA; by default DA=0.25",
    )
    cli.add_argument(
        "-G",
        "--gap-alert",
        metavar="GA",
        type=float,
        default=0.05,
        help="issue an alert in the final report for each marker whose gap rate (proportion of reads containing one or more gap alleles) exceeds DA; by default DA=0.05",
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
        "--hg38",
        default=files("microhapulator") / "data" / "hg38.fasta.gz",
        help=SUPPRESS,
        # Hidden option for testing purposes
    )
    cli.add_argument(
        "--hg38idx",
        default=files("microhapulator") / "data" / "hg38.mmi",
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
