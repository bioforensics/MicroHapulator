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


from hashlib import sha1
from pathlib import Path
from pkg_resources import resource_filename
from subprocess import run
import sys
from tqdm import tqdm
from urllib.request import urlretrieve


class ProgressBar(tqdm):
    """Stolen shamelessly from https://stackoverflow.com/a/53877507/459780."""

    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def compute_shasum(filename):
    with open(filename, "rb") as fh:
        shahash = sha1()
        for block in iter(lambda: fh.read(4096), b""):
            shahash.update(block)
        return shahash.hexdigest()


def download_is_needed(url, filepath, checksum):
    if filepath.exists():
        if not filepath.is_file():
            raise FileNotFoundError(f"{str(filepath)} does not exist or is not a normal file")
        shasum = compute_shasum(filepath)
        if shasum != checksum:
            print(
                f"[MicroHapulator] file {str(filepath)} already present, but checksum failed; re-downloading",
                file=sys.stderr,
            )
        return shasum != checksum
    else:
        return True


def subparser(subparsers):
    desc = "Download and index a GRCh38 assembly file suitable as a whole-genome mapping reference"
    cli = subparsers.add_parser("getrefr", description=desc)


def main(args):
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    checksum = "70fb7af4dff26bffdf27dbef80caf1f0359d488f"
    hg38path = Path(resource_filename("microhapulator", "data/hg38.fasta.gz"))
    if download_is_needed(url, hg38path, checksum):
        print("[MicroHapulator] downloading GRCh38 reference genome", file=sys.stderr)
        with ProgressBar(unit="B", unit_scale=True, miniters=1, desc=hg38path.name) as pb:
            urlretrieve(url, hg38path, reporthook=pb.update_to)
    else:
        print(
            f"[MicroHapulator] file {str(hg38path)} already present, skipping download",
            file=sys.stderr,
        )
    if compute_shasum(hg38path) != checksum:
        raise ValueError(f"checksum failed for {str(hg38path)}")
    index_present = True
    for suffix in ("amb", "ann", "bwt", "pac", "sa"):
        idxfile = Path(f"{str(hg38path)}.{suffix}")
        if not idxfile.is_file():
            index_present = False
            break
    if index_present:
        print(f"[MicroHapulator] BWA index present, good to go!", file=sys.stderr)
    else:
        print(f"[MicroHapulator] building BWA index, this will take some time...", file=sys.stderr)
        run(["bwa", "index", str(hg38path)], check=True)
