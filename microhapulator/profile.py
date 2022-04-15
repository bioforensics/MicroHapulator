# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapulator (https://github.com/bioforensics/microhapulator) and is
# licensed under the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


from collections import defaultdict
from happer.mutate import mutate
from io import StringIO
import json
import jsonschema
from microhapulator import __version__
from microhapulator.parsers import open as mhopen
from microhapulator.parsers import package_file
from numpy.random import choice
import pandas as pd
from pathlib import Path
import sys

SCHEMA = None


def load_schema():
    with mhopen(package_file("data/profile-schema.json"), "r") as fh:
        return json.load(fh)


def check_filter_config(config):
    if config is None:
        return
    expected = set(["Marker", "Static", "Dynamic"])
    missing = expected - set(config)
    if len(missing) > 0:
        missingstr = ",".join(sorted(missing))
        raise ValueError(f"filter config file missing column(s): {missingstr}")
    if len(config.Marker) != len(config.Marker.unique()):
        raise ValueError("filter config file contains duplicate entries for some markers")


class Profile(object):
    def __init__(self, fromfile=None):
        global SCHEMA
        if fromfile:
            if isinstance(fromfile, str) or isinstance(fromfile, Path):
                with mhopen(str(fromfile), "r") as fh:
                    self.data = json.load(fh)
            else:
                self.data = json.load(fromfile)
            if SCHEMA is None:
                SCHEMA = load_schema()
            jsonschema.validate(instance=self.data, schema=SCHEMA)
        else:
            self.data = self.initialize()

    @property
    def gttype(self):
        return "base"

    @property
    def ploidy(self):
        return self.data["ploidy"]

    def initialize(self):
        return {
            "version": __version__,
            "type": self.gttype,
            "ploidy": None,
            "markers": dict(),
        }

    def haploindexes(self):
        hapids = set()
        for markerid, markerdata in self.data["markers"].items():
            for alleledata in markerdata["genotype"]:
                if "index" in alleledata and alleledata["index"] is not None:
                    hapids.add(alleledata["index"])
        assert sorted(hapids) == sorted(range(len(hapids)))
        return hapids

    def markers(self):
        return set(self.data["markers"])

    def haplotypes(self, markerid, index=None):
        if markerid not in self.data["markers"]:
            return set()
        if index is not None:
            return set(
                [
                    a["haplotype"]
                    for a in self.data["markers"][markerid]["genotype"]
                    if "index" in a and a["index"] == index
                ]
            )
        return set([a["haplotype"] for a in self.data["markers"][markerid]["genotype"]])

    def rand_match_prob(self, freqs):
        """Compute the random match probability of this profile.

        Given a set of population allele frequencies, the random match
        probability is the product of the allele frequencies of each allele
        observed in the profile.
        """
        prob = 1.0
        for marker in sorted(self.markers()):
            haplotypes = self.haplotypes(marker)
            diploid_consistent = 1 <= len(haplotypes) <= 2
            if not diploid_consistent:
                nhaps = len(haplotypes)
                msg = f"cannot compute random match prob. for marker with {nhaps} haplotypes"
                raise RandomMatchError(msg)
            result = freqs[(freqs.Marker == marker) & (freqs.Haplotype.isin(haplotypes))]
            if len(haplotypes) == 1:
                p = 0.001
                if len(result) == 1:
                    pp = list(result.Frequency)[0]
                    if pp > 0.0:
                        p = pp
                prob *= p * p
            else:
                p, q = 0.001, 0.001
                if len(result) == 2:
                    pp, qp = list(result.Frequency)
                elif len(result) == 1:
                    pp = list(result.Frequency)[0]
                if pp > 0.0:
                    p = pp
                if len(result) == 2 and qp > 0.0:
                    q = qp
                prob *= 2 * p * q
        return prob

    def rmp_lr_test(self, other, freqs, erate=0.001):
        """Compute a likelihood ratio test for the random match probability.

        The likelihood ratio test compares the probability that the two samples
        come from the same source to the probability that the samples come from
        two unrelated sources. The first probability (numerator) should be 1.0,
        with any missing or incongruent alleles treated as genotyping error.
        The second probability (denominator) is the random match probability.
        The test only makes sense when the two profiles being compared are
        identical or nearly identical.
        """
        assert self.data["ploidy"] in (2, None)
        assert other.data["ploidy"] in (2, None)
        mismatches = 0
        for marker in self.markers():
            selfhaps = self.haplotypes(marker)
            otherhaps = other.haplotypes(marker)
            if selfhaps == otherhaps:
                pass
            elif len(selfhaps & otherhaps) == 1:
                mismatches += 1
            else:
                mismatches += 2
        numerator = erate**mismatches
        denominator = self.rand_match_prob(freqs)
        return numerator / denominator

    def dump(self, outfile):
        if isinstance(outfile, str) or isinstance(outfile, Path):
            with mhopen(str(outfile), "w") as fh:
                json.dump(self.data, fh, indent=4, sort_keys=True)
        else:
            json.dump(self.data, outfile, indent=4, sort_keys=True)

    def __eq__(self, other):
        if not isinstance(other, Profile):
            return False
        if self.markers() != other.markers():
            return False
        for marker in self.markers():
            if self.haplotypes(marker) != other.haplotypes(marker):
                return False
        return True

    def __str__(self):
        return json.dumps(self.data, indent=4, sort_keys=True)

    def bedstream(self, markers):
        for marker in sorted(self.markers()):
            result = markers[markers.Marker == marker]
            if len(result) == 0:
                raise ValueError(f"unknown marker identifier '{marker}'")
            offsets = sorted(result.Offset)
            variants = [list() for _ in range(len(offsets))]
            for index in sorted(self.haploindexes()):
                haplotype = self.haplotypes(marker, index=index).pop()
                for snp, allelelist in zip(haplotype.split(","), variants):
                    allelelist.append(snp)
            for offset, snps in zip(offsets, variants):
                haplostr = "|".join(snps)
                yield "\t".join((marker, str(offset), str(offset + 1), haplostr))

    def seqstream(self, refrseqs):
        for marker in sorted(self.markers()):
            yield marker, refrseqs[marker]

    def bedstr(self, markers):
        out = StringIO()
        for line in self.bedstream(markers):
            print(line, file=out)
        return out.getvalue()

    def haploseqs(self, markers, refrseqs):
        """Apply genotype to reference and construct full haplotype sequences."""
        ss = self.seqstream(refrseqs)
        bs = self.bedstream(markers)
        for defline, sequence in mutate(ss, bs):
            yield defline, sequence

    def unite(mom, dad):
        """Simulate the creation of a new profile from a mother and father

        :param microhapulator.profile.Profile mom: typing result or simulated profile
        :param microhapulator.profile.Profile dad: typing result or simulated profile
        :returns: a simulated offspring
        :rtype: microhapulator.profile.SimulatedProfile
        """
        prof = SimulatedProfile(ploidy=2)
        allmarkers = mom.markers() | dad.markers()
        commonmarkers = mom.markers() & dad.markers()
        if len(commonmarkers) == 0:
            raise ValueError("mom and dad profiles have no markers in common")
        notshared = allmarkers - commonmarkers
        if len(notshared) > 0:
            message = "markers not common to mom and dad profiles are excluded: "
            message += ", ".join(notshared)
            print("[MicroHapulator::profile]", message, file=sys.stderr)
        for parent, hapid in zip((mom, dad), (0, 1)):
            for marker in sorted(commonmarkers):
                haplotype = choice(sorted(parent.haplotypes(marker)))
                prof.add(hapid, marker, haplotype)
        return prof

    def unmix(self):
        assert self.ploidy % 2 == 0
        ncontrib = int(self.ploidy / 2)
        profiles = [SimulatedProfile() for _ in range(ncontrib)]
        for marker in self.markers():
            for contrib in range(ncontrib):
                for subindex in range(2):
                    index = (2 * contrib) + subindex
                    haplotype = self.haplotypes(marker, index=index).pop()
                    profiles[contrib].add(subindex, marker, haplotype)
        return profiles


class SimulatedProfile(Profile):
    """Profile represented by phased alleles at a number of specified markers.

    Alleles are stored in a dictionary, with microhap marker ID/name as the key
    and a list as the value. Typically each list contains 2 items, the
    haplotype/phase 0 allele and the hap/phase 1 allele. However, if the
    profile represents a mixture there may be more than 2 haplotypes, as
    indicated by the `ploidy` parameter.
    """

    def populate_from_bed(bedfile):
        with mhopen(bedfile, "r") as fh:
            line = next(fh)
            ploidy = line.count("|") + 1
            fh.seek(0)
            marker_alleles = defaultdict(lambda: [list() for _ in range(ploidy)])
            for line in fh:
                line = line.strip()
                if line == "":
                    continue
                marker, start, end, allelestr = line.split("\t")
                alleles = allelestr.split("|")
                for i, a in enumerate(alleles):
                    marker_alleles[marker][i].append(a)
            profile = SimulatedProfile(ploidy=ploidy)
            for marker, allele_list in marker_alleles.items():
                for i, haplotype in enumerate(allele_list):
                    profile.add(i, marker, ",".join(haplotype))
            return profile

    def merge(profiles):
        """Combine simulated profiles into a mock DNA mixture

        :param list profiles: list of simulated profiles
        :returns: a combined profile
        :rtype: microhapulator.profile.SimulatedProfile
        """
        ploidy = 2 * len(profiles)
        prof = SimulatedProfile(ploidy=ploidy)
        offset = 0
        for profile in profiles:
            for markerid, markerdata in sorted(profile.data["markers"].items()):
                for haplotype in markerdata["genotype"]:
                    prof.add(offset + haplotype["index"], markerid, haplotype["haplotype"])
            offset += 2
        return prof

    def __init__(self, fromfile=None, ploidy=0):
        super(SimulatedProfile, self).__init__(fromfile=fromfile)
        if ploidy:
            self.data["ploidy"] = ploidy

    def add(self, hapid, marker, haplotype):
        if self.data["ploidy"] is not None and self.data["ploidy"] > 0:
            assert hapid in range(self.data["ploidy"])
        if marker not in self.data["markers"]:
            self.data["markers"][marker] = {"genotype": list()}
        self.data["markers"][marker]["genotype"].append({"haplotype": haplotype, "index": hapid})

    @property
    def gttype(self):
        return "SimulatedProfile"


class TypingResult(Profile):
    def __init__(self, fromfile=None):
        super(TypingResult, self).__init__(fromfile=fromfile)

    def record_coverage(self, marker, cov_by_pos, ndiscarded=0):
        self.data["markers"][marker] = {
            "genotype": list(),
            "mean_coverage": 0.0,
            "min_coverage": 0,
            "max_coverage": 0,
            "num_discarded_reads": ndiscarded,
            "typing_result": dict(),
        }
        if len(cov_by_pos) > 0:
            avgcov = sum(cov_by_pos) / len(cov_by_pos)
            self.data["markers"][marker]["mean_coverage"] = round(avgcov, 1)
            self.data["markers"][marker]["min_coverage"] = min(cov_by_pos)
            self.data["markers"][marker]["max_coverage"] = max(cov_by_pos)

    def record_haplotype(self, marker, haplotype, count):
        self.data["markers"][marker]["typing_result"][haplotype] = count

    def filter(self, static=None, dynamic=None, config=None):
        """Apply static and/or dynamic thresholds to distinguish true and false haplotypes

        Thresholds are applied to the haplotype read counts of a raw typing result. Static integer
        thresholds are commonly used as detection thresholds, below which any haplotype count is
        considered noise. Dynamic thresholds are commonly used as analytical thresholds and
        represent a percentage of the total read count at the marker, after any haplotypes failing
        a static threshold are discarded.

        :param int static: global fixed read count threshold
        :param float dynamic: global percentage of total read count; e.g. use `dynamic=0.02` to apply a 2% analytical threshold
        :param pandas.DataFrame config: tabular data structure specifying marker-specific thresholds to override global thresholds; three required columns: **Marker** for the marker name; **Static** and **Dynamic** for marker-specific thresholds
        """
        check_filter_config(config)
        for marker, mdata in self.data["markers"].items():
            markerstatic = static
            markerdynamic = dynamic
            if config is not None:
                thresh = config[config.Marker == marker]
                assert len(thresh) in (0, 1)
                if len(thresh) == 1:
                    markerstatic = thresh.Static.iloc[0]
                    markerdynamic = thresh.Dynamic.iloc[0]
            self.data["markers"][marker]["genotype"] = list()
            if markerstatic is None and markerdynamic is None:
                # No thresholds for calling haplotypes, just report raw haplotype counts
                continue
            hapcounts = mdata["typing_result"]
            genotype_call = set()
            filtered = set()
            totalcount = sum(hapcounts.values())
            filteredcount = 0
            if markerstatic is not None and markerstatic > 0:
                for haplotype, count in hapcounts.items():
                    if count < markerstatic:
                        filtered.add(haplotype)
                        filteredcount += count
            if markerdynamic is not None and markerdynamic > 0.0:
                for haplotype, count in hapcounts.items():
                    if haplotype in filtered:
                        continue
                    if count < (totalcount - filteredcount) * markerdynamic:
                        filtered.add(haplotype)
                    else:
                        genotype_call.add(haplotype)
            self.data["markers"][marker]["genotype"] = [
                {"haplotype": ht} for ht in sorted(genotype_call)
            ]

    def dump_csv(self, outfile, samplename, counts=True, fix_homo=False):
        max_haps = 0
        for marker in sorted(self.markers()):
            max_haps = max(max_haps, len(self.haplotypes(marker)))
        entries = list()
        column_names = ["SampleName", "Marker"]
        column_names += [f"Allele{i+1}" for i in range(max_haps)]
        if counts is True:
            column_names += [f"Height{i+1}" for i in range(max_haps)]
        for marker in sorted(self.markers()):
            haplotypes = sorted(self.haplotypes(marker))
            if counts is False and len(haplotypes) == 1 and fix_homo:
                haplotypes = haplotypes * 2
            hapcounts = list()
            if counts is True:
                for ht in haplotypes:
                    count = self.data["markers"][marker]["typing_result"][ht]
                    hapcounts.append(count)
                while len(hapcounts) < max_haps:
                    hapcounts.append(None)
            # May need to convert , to - pending a test of LRMix and EuroForMix's ability to read
            # CSV files with quoted strings containing "," characters.
            # haplotypes = [ht.replace(",", "-") for ht in haplotypes]
            while len(haplotypes) < max_haps:
                haplotypes.append(None)
            entry = [samplename, marker, *haplotypes, *hapcounts]
            entries.append(entry)
        table = pd.DataFrame(entries, columns=column_names)
        table.to_csv(outfile, index=False)

    def typing_rate(self):
        data = {
            "Marker": list(),
            "TypedReads": list(),
            "TotalReads": list(),
            "TypingRate": list(),
        }
        for marker, mdata in self.data["markers"].items():
            num_typed_reads = 0
            for mhallele, count in mdata["typing_result"].items():
                num_typed_reads += count
            total_reads = num_typed_reads + mdata["num_discarded_reads"]
            rate = 0.0
            if total_reads > 0:
                rate = num_typed_reads / total_reads
            data["Marker"].append(marker)
            data["TypedReads"].append(num_typed_reads)
            data["TotalReads"].append(total_reads)
            data["TypingRate"].append(rate)
        return pd.DataFrame(data)

    @property
    def gttype(self):
        return "TypingResult"
