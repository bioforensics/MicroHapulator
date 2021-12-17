#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapulator (github.com/bioforensics/microhapulator)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
from happer.mutate import mutate
from io import StringIO
import json
import jsonschema
import microhapdb
from microhapdb.marker import TargetAmplicon
import microhapulator
from numpy.random import choice


def load_schema():
    with microhapulator.open(microhapulator.package_file("profile-schema.json"), "r") as fh:
        return json.load(fh)


schema = None


class Profile(object):
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
            microhapulator.plog("[MicroHapulator::profile]", message)
        for parent, hapid in zip((mom, dad), (0, 1)):
            for marker in sorted(commonmarkers):
                haploallele = choice(sorted(parent.alleles(marker)))
                gt.add(hapid, marker, haploallele)
        return gt

    def __init__(self, fromfile=None):
        global schema
        if fromfile:
            if isinstance(fromfile, str):
                with microhapulator.open(fromfile, "r") as fh:
                    self.data = json.load(fh)
            else:
                self.data = json.load(fromfile)
            if schema is None:
                schema = load_schema()
            jsonschema.validate(instance=self.data, schema=schema)
        else:
            self.data = self.initialize()

    @property
    def gttype(self):
        return "base"

    @property
    def ploidy(self):
        return self.data["ploidy"]

    def initialize(self):
        data = {
            "version": microhapulator.__version__,
            "type": self.gttype,
            "ploidy": None,
            "markers": dict(),
        }
        return data

    def haplotypes(self):
        hapids = set()
        for markerid, markerdata in self.data["markers"].items():
            for alleledata in markerdata["genotype"]:
                if "haplotype" in alleledata:
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
            result = freqs[(freqs.Marker == marker) & (freqs.Allele.isin(alleles))]
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
        if isinstance(outfile, str):
            with microhapulator.open(outfile, "w") as fh:
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

    @property
    def bedstream(self):
        hapids = self.haplotypes()
        for marker in sorted(self.markers()):
            result = microhapdb.markers[microhapdb.markers.Name == marker]
            if len(result) == 0:
                raise ValueError('unknown marker identifier "{:s}"'.format(marker))
            markerdata = result.iloc[0]
            context = TargetAmplicon(markerdata, delta=30, minlen=350)
            coords = list(map(int, markerdata.Offsets.split(",")))
            coords = list(map(context.global_to_local, coords))
            variants = [list() for _ in range(len(coords))]
            for haplotype in sorted(hapids):
                allele = self.alleles(marker, haplotype=haplotype).pop()
                for var, varlist in zip(allele.split(","), variants):
                    varlist.append(var)
            for coord, var in zip(coords, variants):
                allelestr = "|".join(var)
                yield "\t".join((marker, str(coord), str(coord + 1), allelestr))

    @property
    def seqstream(self):
        for marker in sorted(self.markers()):
            canonid = microhapdb.markers[microhapdb.markers.Name == marker].iloc[0]
            amp = TargetAmplicon(canonid, delta=30, minlen=350)
            yield amp.defline, amp.amplicon_seq

    @property
    def bedstr(self):
        out = StringIO()
        for line in self.bedstream:
            print(line, file=out)
        return out.getvalue()

    @property
    def haploseqs(self):
        """Apply genotype to reference and construct full haplotype sequences."""
        mutator = mutate(self.seqstream, self.bedstream)
        for defline, sequence in mutator:
            yield defline, sequence

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

    >>> gt = SimulatedProfile()
    >>> gt.add(0, 'mh21KK-315', 'G,C,T')
    >>> gt.add(1, 'mh21KK-315', 'A,T,C')
    >>> gt.add(0, 'mh21KK-316', 'A,C,G,T')
    >>> gt.add(1, 'mh21KK-316', 'A,T,G,C')
    >>> print(gt.bedstr)
    mh21KK-315 102     103     G|A
    mh21KK-315 207     208     C|T
    mh21KK-315 247     248     T|C
    mh21KK-316 108     109     A|A
    mh21KK-316 132     133     C|T
    mh21KK-316 179     180     G|G
    mh21KK-316 242     243     T|C
    """

    def populate_from_bed(bedfile):
        with microhapulator.open(bedfile, "r") as fh:
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


class ObservedProfile(Profile):
    def __init__(self, fromfile=None):
        super(ObservedProfile, self).__init__(fromfile=fromfile)

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
        return "ObservedProfile"
