# find-ibs-segments
#
# Copyright (C) 2026 LaKisha T. David
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Tests of the run finder and of the PBWT it uses.

Run from the repository root: ``python -m unittest discover -s tests``.

The PBWT is compared with a direct comparison of every pair of haplotypes on
seeded random data.  The hand-built cases use markers 1 kb and 0.125 cM
apart (an exact binary fraction, so lengths are exact) and state their
expected output.  All data are synthetic and built here.
"""

from __future__ import annotations

import os
import random
import subprocess
import sys
import tempfile
import unittest
from typing import Dict, List, Sequence, Tuple

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from find_ibs_segments import find  # noqa: E402
from find_ibs_segments.java_compat import JavaTermination  # noqa: E402
from find_ibs_segments.pbwt import long_matches  # noqa: E402

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
POSITIONS = [1000 * i for i in range(1, 41)]
STEP_CM = 0.125


def write_inputs(directory: str, genotypes: Dict[Tuple[str, int], Dict[str, str]],
                 samples: Sequence[str] = ("A", "B"), chroms: Sequence[str] = ("1",),
                 positions: Sequence[int] = POSITIONS,
                 records: Sequence[str] = (),
                 alts: Dict[Tuple[str, int], str] = None) -> Tuple[str, str]:
    """A phased VCF (0|0 unless ``genotypes[(chrom, pos)][sample]`` says
    otherwise; ALT G unless ``alts`` says otherwise) and a map with a point
    at every marker."""
    lines = ["##fileformat=VCFv4.2", HEADER + "".join("\t" + s for s in samples)]
    for c in chroms:
        for p in positions:
            gts = [genotypes.get((c, p), {}).get(s, "0|0") for s in samples]
            alt = (alts or {}).get((c, p), "G")
            lines.append("\t".join([c, str(p), ".", "A", alt, ".", "PASS", ".",
                                    "GT"] + gts))
    lines.extend(records)
    vcf = os.path.join(directory, "in.vcf")
    with open(vcf, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    mp = os.path.join(directory, "in.map")
    with open(mp, "w", encoding="utf-8", newline="\n") as f:
        for c in chroms:
            for i, p in enumerate(sorted(positions)):
                f.write("%s\t.\t%r\t%d\n" % (c, STEP_CM * i, p))
    return vcf, mp


def run(genotypes, *extra, **kw) -> find.Result:
    with tempfile.TemporaryDirectory() as d:
        vcf, mp = write_inputs(d, genotypes, **kw)
        return find.run([vcf, mp] + list(extra))


def stopped(genotypes, *extra, **kw) -> JavaTermination:
    try:
        run(genotypes, *extra, **kw)
    except JavaTermination as t:
        return t
    raise AssertionError("the program did not stop")


def row(s1, h1, s2, h2, start, end, cm, markers, carriers, freq, chrom="1"):
    return "\t".join(str(x) for x in (s1, h1, s2, h2, chrom, start, end, cm,
                                      markers, carriers, freq))


# A's haplotype 1 carries 1 at 21000; B's haplotype 2 carries 1 at 35000.
GENO = {("1", 21000): {"A": "1|0"}, ("1", 35000): {"B": "0|1"}}
EXPECTED = [
    row("A", 1, "A", 2, 1000, 20000, "2.375", 20, 4, "1"),
    row("A", 1, "B", 1, 1000, 20000, "2.375", 20, 4, "1"),
    row("A", 1, "B", 2, 1000, 20000, "2.375", 20, 4, "1"),
    row("A", 2, "B", 2, 1000, 34000, "4.125", 34, 3, "0.75"),
    row("B", 1, "B", 2, 1000, 34000, "4.125", 34, 3, "0.75"),
    row("A", 2, "B", 1, 1000, 40000, "4.875", 40, 2, "0.5"),
    row("A", 1, "B", 2, 22000, 34000, "1.5", 13, 4, "1"),
    row("A", 1, "A", 2, 22000, 40000, "2.25", 19, 3, "0.75"),
    row("A", 1, "B", 1, 22000, 40000, "2.25", 19, 3, "0.75"),
]


class Runs(unittest.TestCase):

    def test_every_pair_and_its_columns(self):
        # Pairs of one sample's two haplotypes are included; a run ends at
        # the marker before a discordant one, or at the chromosome's end.
        got = run(GENO, "1.5")
        self.assertEqual(got.output.splitlines(), EXPECTED)
        self.assertEqual(got.summary,
                         "find-ibs-segments: samples 2, markers 40, "
                         "chromosomes 1; runs of at least 1.5 cM: 9 (3 between "
                         "a sample's own two haplotypes)\n")

    def test_minimum_length_boundary_and_default(self):
        # A1-B2 22000-34000 is exactly 1.5 cM; 1.5000000000000002 is the
        # next double above 1.5.
        got = run(GENO, "1.5000000000000002").output.splitlines()
        self.assertEqual(got, [r for r in EXPECTED if "\t22000\t34000\t" not in r])
        # Without [min_cm] the minimum is 3 cM.
        self.assertEqual(run(GENO).output.splitlines(),
                         [r for r in EXPECTED
                          if float(r.split("\t")[7]) >= 3.0])

    def test_runs_stay_on_their_chromosome(self):
        geno = {("2", 1000): {"A": "1|0"}}
        got = run(geno, "4", chroms=("1", "2")).output.splitlines()
        whole = [row(a, h, b, k, 1000, 40000, "4.875", 40, 4, "1", chrom="1")
                 for a, h, b, k in (("A", 1, "A", 2), ("A", 1, "B", 1),
                                    ("A", 1, "B", 2), ("A", 2, "B", 1),
                                    ("A", 2, "B", 2), ("B", 1, "B", 2))]
        # On chromosome 2, A's haplotype 1 differs at the first marker.
        second = [row(a, h, b, k, 2000, 40000, "4.75", 39, 4, "1", chrom="2")
                  for a, h, b, k in (("A", 1, "A", 2), ("A", 1, "B", 1),
                                     ("A", 1, "B", 2))]
        second += [row(a, h, b, k, 1000, 40000, "4.875", 40, 3, "0.75", chrom="2")
                   for a, h, b, k in (("A", 2, "B", 1), ("A", 2, "B", 2),
                                      ("B", 1, "B", 2))]
        self.assertEqual(got, whole + sorted(second, key=lambda r: (
            int(r.split("\t")[5]), r.split("\t")[0], r.split("\t")[1],
            r.split("\t")[2], r.split("\t")[3])))

    def test_frequency_format(self):
        for carriers, total, text in ((2, 96, "0.0208333"), (1, 3, "0.333333"),
                                      (8, 8, "1"), (1, 2, "0.5"),
                                      (2, 30000, "0.0000666667")):
            with self.subTest(carriers=carriers, total=total):
                self.assertEqual(find.format_frequency(carriers, total), text)


class Errors(unittest.TestCase):

    def test_biallelic_markers_only(self):
        t = stopped({("1", 5000): {"A": "2|0"}}, "1.5", alts={("1", 5000): "G,T"})
        self.assertEqual(t.stderr, "ERROR: find-ibs-segments needs biallelic "
                                   "markers: chromosome 1 position 5000 has "
                                   "allele index 2\n")

    def test_records_in_position_order_in_one_block_per_chromosome(self):
        swapped = POSITIONS[:]
        swapped[4], swapped[5] = swapped[5], swapped[4]
        t = stopped({}, "1.5", positions=swapped)
        self.assertEqual(t.stderr, "ERROR: VCF records are not in position "
                                   "order: chromosome 1 position 5000 follows "
                                   "position 6000\n")
        extra = ["1\t50000\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0"]
        t = stopped({}, "1.5", chroms=("1", "2"), records=extra)
        self.assertEqual(t.stderr, "ERROR: the records of chromosome 1 are not "
                                   "in one block in the VCF file\n")

    def test_arguments(self):
        with self.assertRaises(JavaTermination) as ctx:
            find.run([])
        self.assertEqual((ctx.exception.status, ctx.exception.stdout),
                         (0, find.usage() + "\n"))
        for args, message in ((["a"], "ERROR: expected 2 or 3 arguments, found 1"),
                              (["a", "b", "1", "x"],
                               "ERROR: expected 2 or 3 arguments, found 4")):
            with self.subTest(args=args):
                with self.assertRaises(JavaTermination) as ctx:
                    find.run(args)
                self.assertEqual(ctx.exception.status, 1)
                self.assertEqual(ctx.exception.stderr,
                                 find.usage() + "\n" + message + "\n")
        for value, message in (("x", "ERROR: [min_cm] is not a numerical value: x"),
                               ("0", "ERROR: [min_cm]=0"), ("-1", "ERROR: [min_cm]=-1"),
                               ("NaN", "ERROR: [min_cm]=NaN")):
            with self.subTest(min_cm=value):
                t = stopped({}, value)
                self.assertEqual(t.stderr, find.usage() + "\n" + message + "\n")

    def test_command_line(self):
        with tempfile.TemporaryDirectory() as d:
            vcf, mp = write_inputs(d, GENO)
            done = subprocess.run([sys.executable, "-m", "find_ibs_segments", vcf, mp,
                                   "1.5"], capture_output=True, cwd=REPO)
        self.assertEqual(done.returncode, 0)
        self.assertEqual(done.stdout.decode().splitlines(), EXPECTED)
        self.assertTrue(done.stderr.decode().startswith("find-ibs-segments: "))


def brute_force(sites: List[bytes], cm: List[float],
                min_cm: float) -> set:
    """Every maximal identical run of every pair, found directly."""
    n = len(sites)
    found = set()
    for x in range(len(sites[0])):
        for y in range(x + 1, len(sites[0])):
            k = 0
            while k < n:
                if sites[k][x] != sites[k][y]:
                    k += 1
                    continue
                start = k
                while k < n and sites[k][x] == sites[k][y]:
                    k += 1
                if cm[k - 1] - cm[start] >= min_cm:
                    found.add((x, y, start, k - 1))
    return found


class PbwtAgainstBruteForce(unittest.TestCase):

    def test_seeded_random_cases(self):
        total = 0
        for seed in range(300):
            with self.subTest(seed=seed):
                rng = random.Random(seed)
                n_haps, n_sites = rng.randint(2, 14), rng.randint(1, 60)
                founders = [[rng.randint(0, 1) for _ in range(n_sites)]
                            for _ in range(rng.randint(1, 4))]
                haps = []
                for _ in range(n_haps):
                    f = rng.randrange(len(founders))
                    h = []
                    for k in range(n_sites):
                        if rng.random() < 0.1:
                            f = rng.randrange(len(founders))
                        allele = founders[f][k]
                        if rng.random() < 0.03:
                            allele = 1 - allele
                        h.append(allele)
                    haps.append(h)
                sites = [bytes(haps[h][k] for h in range(n_haps))
                         for k in range(n_sites)]
                cm, c = [], 0.0
                for _ in range(n_sites):
                    c += 0.0 if rng.random() < 0.2 else rng.uniform(0, 0.5)
                    cm.append(c)
                min_cm = rng.choice([0.01, 0.3, 1.0, 2.5])
                got = list(long_matches(sites, cm, min_cm))
                self.assertEqual(len(got), len(set(got)))
                self.assertEqual(set(got), brute_force(sites, cm, min_cm))
                total += len(got)
        self.assertGreater(total, 1000)


if __name__ == "__main__":
    unittest.main()
