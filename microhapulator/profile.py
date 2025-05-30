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
from importlib.resources import files
from io import StringIO
import json
import jsonschema
from microhapulator import __version__
from microhapulator import open as mhopen
from microhapulator.happer.mutate import mutate
from numpy.random import choice
import pandas as pd
from pathlib import Path
import sys

SCHEMA = None


class RandomMatchError(ValueError):
    pass


def load_schema():
    with mhopen(files("microhapulator") / "data/profile-schema.json", "r") as fh:
        return json.load(fh)


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

    def loci(self):
        return set(m.split(".")[0] for m in self.markers())

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

    def rand_match_prob(self, freqs, minfreq=0.01):
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
            result = freqs[(freqs.Marker == marker) & (freqs.Allele.isin(haplotypes))]
            if len(haplotypes) == 1:
                p = minfreq
                if len(result) == 1:
                    pp = list(result.Frequency)[0]
                    if pp > 0.0:
                        p = pp
                prob *= p * p
            else:
                p, q = minfreq, minfreq
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

    def bedstream(self, mhindex):
        mhindex.validate(refrids=self.loci(), symmetric=True)
        for markerid in sorted(self.markers()):
            marker = mhindex.markers[markerid]
            offsets = marker.offsets_locus
            variants = [list() for _ in range(len(offsets))]
            for i in sorted(self.haploindexes()):
                haplotype = self.haplotypes(markerid, index=i).pop()
                for snp, allelelist in zip(haplotype.split(":"), variants):
                    allelelist.append(snp)
            for offset, snps in zip(offsets, variants):
                haplostr = "|".join(snps)
                yield "\t".join((markerid, str(offset), str(offset + 1), haplostr))

    def seqstream(self, mhindex):
        for markerid in sorted(self.markers()):
            yield markerid, mhindex.sequence(markerid)

    def bedstr(self, mhindex):
        out = StringIO()
        for line in self.bedstream(mhindex):
            print(line, file=out)
        return out.getvalue()

    def haploseqs(self, mhindex):
        """Apply genotype to reference and construct full haplotype sequences."""
        ss = self.seqstream(mhindex)
        bs = self.bedstream(mhindex)
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
                    profile.add(i, marker, ":".join(haplotype))
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

    def filter(self, thresholds):
        """Apply static and/or dynamic thresholds to distinguish true and false haplotypes

        Thresholds are applied to the haplotype read counts of a raw typing result. Static integer
        thresholds are commonly used as detection thresholds, below which any haplotype count is
        considered noise. Dynamic thresholds are commonly used as analytical thresholds and
        represent a percentage of the total read count at the marker, after any haplotypes failing
        a static threshold are discarded.

        :param ThresholdIndex thresholds: data structure containing static and dynamic thresholds
        """
        for marker, mdata in self.data["markers"].items():
            self.data["markers"][marker]["genotype"] = list()
            self.data["markers"][marker]["thresholds"] = dict()
            hapcounts = mdata["typing_result"]
            genotype_call = set()
            filtered = set()
            totalcount = sum(hapcounts.values())
            filteredcount = 0
            static, dynamic = thresholds.get(marker)
            if static is not None and not pd.isna(static) and static > 0:
                self.data["markers"][marker]["thresholds"]["static"] = int(static)
                for haplotype, count in hapcounts.items():
                    if count < static:
                        filtered.add(haplotype)
                        filteredcount += count
            if dynamic is not None and not pd.isna(dynamic) and dynamic > 0.0:
                threshold = (totalcount - filteredcount) * dynamic
                self.data["markers"][marker]["thresholds"]["dynamic"] = threshold
                for haplotype, count in hapcounts.items():
                    if count < threshold:
                        filtered.add(haplotype)
            for haplotype, count in hapcounts.items():
                if haplotype not in filtered:
                    genotype_call.add(haplotype)
            self.data["markers"][marker]["genotype"] = [
                {"haplotype": ht} for ht in sorted(genotype_call)
            ]

    def allele_counts(self, marker):
        if marker not in self.data["markers"]:
            raise IndexError(marker)
        return self.data["markers"][marker]["typing_result"]

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
        table.to_csv(outfile, index=False, float_format="%d")

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

    def gap_rate(self):
        data = []
        for marker, mdata in self.data["markers"].items():
            num_typed_reads, num_gap_reads = 0, 0
            for mhallele, count in mdata["typing_result"].items():
                num_typed_reads += count
                if "-" in mhallele:
                    num_gap_reads += count
            rate = 0.0
            if num_typed_reads > 0 and num_gap_reads > 0:
                rate = num_gap_reads / num_typed_reads
            entry = (marker, num_typed_reads, num_gap_reads, rate)
            data.append(entry)
        return pd.DataFrame(data, columns=["Marker", "TypedReads", "GappedReads", "GappedRate"])

    @property
    def gttype(self):
        return "TypingResult"
