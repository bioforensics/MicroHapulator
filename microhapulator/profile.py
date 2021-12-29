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
from pathlib import Path
import sys

SCHEMA = None


def load_schema():
    with mhopen(package_file("profile-schema.json"), "r") as fh:
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

    def haplotypes(self):
        hapids = set()
        for markerid, markerdata in self.data["markers"].items():
            for alleledata in markerdata["genotype"]:
                if "haplotype" in alleledata and alleledata["haplotype"] is not None:
                    hapids.add(alleledata["haplotype"])
        assert sorted(hapids) == sorted(range(len(hapids)))
        return hapids

    def markers(self):
        return set(list(self.data["markers"]))

    def alleles(self, markerid, haplotype=None):
        if markerid not in self.data["markers"]:
            return set()
        if haplotype is not None:
            return set(
                [
                    a["allele"]
                    for a in self.data["markers"][markerid]["genotype"]
                    if "haplotype" in a and a["haplotype"] == haplotype
                ]
            )
        return set([a["allele"] for a in self.data["markers"][markerid]["genotype"]])

    def rand_match_prob(self, freqs):
        """Compute the random match probability of this profile.

        Given a set of population allele frequencies, the random match
        probability is the product of the allele frequencies of each allele
        observed in the profile.
        """
        prob = 1.0
        for marker in sorted(self.markers()):
            alleles = self.alleles(marker)
            diploid_consistent = 1 <= len(alleles) <= 2
            if not diploid_consistent:
                msg = f"cannot compute random match prob. for marker with {len(alleles)} alleles"
                raise RandomMatchError(msg)
            result = freqs[(freqs.Marker == marker) & (freqs.Haplotype.isin(alleles))]
            if len(alleles) == 1:
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
            selfalleles = self.alleles(marker)
            otheralleles = other.alleles(marker)
            if selfalleles == otheralleles:
                pass
            elif len(selfalleles & otheralleles) == 1:
                mismatches += 1
            else:
                mismatches += 2
        numerator = erate ** mismatches
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
            if self.alleles(marker) != other.alleles(marker):
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
            for haplotype in sorted(self.haplotypes()):
                allele = self.alleles(marker, haplotype=haplotype).pop()
                for var, varlist in zip(allele.split(","), variants):
                    varlist.append(var)
            for offset, var in zip(offsets, variants):
                haplostr = "|".join(var)
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
        """Simulate the creation of a new profile from a mother and father."""
        gt = SimulatedProfile(ploidy=2)
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
                haploallele = choice(sorted(parent.alleles(marker)))
                gt.add(hapid, marker, haploallele)
        return gt

    def unmix(self):
        assert self.ploidy % 2 == 0
        ncontrib = int(self.ploidy / 2)
        profiles = [SimulatedProfile() for _ in range(ncontrib)]
        for marker in self.markers():
            for contrib in range(ncontrib):
                for hap in range(2):
                    haplotype = (2 * contrib) + hap
                    allele = self.alleles(marker, haplotype=haplotype).pop()
                    profiles[contrib].add(hap, marker, allele)
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
                for i, allele in enumerate(allele_list):
                    profile.add(i, marker, ",".join(allele))
            return profile

    def merge(profiles):
        ploidy = 2 * len(profiles)
        gt = SimulatedProfile(ploidy=ploidy)
        offset = 0
        for profile in profiles:
            for markerid, markerdata in sorted(profile.data["markers"].items()):
                for allele in markerdata["genotype"]:
                    gt.add(offset + allele["haplotype"], markerid, allele["allele"])
            offset += 2
        return gt

    def __init__(self, fromfile=None, ploidy=0):
        super(SimulatedProfile, self).__init__(fromfile=fromfile)
        if ploidy:
            self.data["ploidy"] = ploidy

    def add(self, hapid, marker, allele):
        if self.data["ploidy"] is not None and self.data["ploidy"] > 0:
            assert hapid in range(self.data["ploidy"])
        if marker not in self.data["markers"]:
            self.data["markers"][marker] = {"genotype": list()}
        self.data["markers"][marker]["genotype"].append({"allele": allele, "haplotype": hapid})

    @property
    def gttype(self):
        return "SimulatedProfile"


class TypingResult(Profile):
    def __init__(self, fromfile=None):
        super(TypingResult, self).__init__(fromfile=fromfile)

    def record_coverage(self, marker, cov_by_pos, ndiscarded=0):
        self.data["markers"][marker] = {
            "mean_coverage": 0.0,
            "min_coverage": 0,
            "max_coverage": 0,
            "num_discarded_reads": ndiscarded,
            "allele_counts": dict(),
        }
        if len(cov_by_pos) > 0:
            avgcov = sum(cov_by_pos) / len(cov_by_pos)
            self.data["markers"][marker]["mean_coverage"] = round(avgcov, 1)
            self.data["markers"][marker]["min_coverage"] = min(cov_by_pos)
            self.data["markers"][marker]["max_coverage"] = max(cov_by_pos)

    def record_allele(self, marker, allele, count):
        self.data["markers"][marker]["allele_counts"][allele] = count

    def infer(self, ecthreshold=0.25, static=None, dynamic=None):
        for marker, mdata in self.data["markers"].items():
            if static is None and dynamic is None:
                # No thresholds for calling haplotypes, just report raw haplotype counts
                self.data["markers"][marker]["genotype"] = list()
                continue
            gt = set()
            allelecounts = mdata["allele_counts"]
            for allele, count in allelecounts.items():
                eff_cov = 1.0 - (mdata["num_discarded_reads"] / mdata["max_coverage"])
                if dynamic is None or (static is not None and eff_cov < ecthreshold):
                    # Use static cutoff (low effective coverage, or dynamic cutoff undefined)
                    if count < static:
                        continue
                else:
                    # Use dynamic cutoff (high effective coverage, or static cutoff undefined)
                    if len(allelecounts.values()) == 0:
                        continue
                    avgcount = sum(allelecounts.values()) / len(allelecounts.values())
                    if count < avgcount * dynamic:
                        continue
                gt.add(allele)
            self.data["markers"][marker]["genotype"] = [
                {"allele": a, "haplotype": None} for a in sorted(gt)
            ]

    @property
    def gttype(self):
        return "TypingResult"
