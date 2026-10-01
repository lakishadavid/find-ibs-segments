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
"""The program: lists every maximal run of identical alleles, at least
[min_cm] cM long, between two haplotypes of a phased VCF file.

Usage: ``python -m find_ibs_segments [vcf] [map] [min_cm] > [out]``, where
``[min_cm]`` is optional (default 3).

Every pair of haplotypes is compared: haplotypes of two samples, and a
sample's own two haplotypes.  A run is a stretch of consecutive markers at
which the two haplotypes carry the same allele, extended in both directions
until the next marker differs or the chromosome's markers end, as Refined
IBD extends its segments (IbsHapSegments.inclStart and inclEnd in
Refined-IBD 17Jan20.102).  Its length is the genetic map distance between its
first and last marker, and it is kept if that is at least [min_cm], as
Refined IBD's length test does.  There is no LOD score and no analysis
window: each chromosome is searched whole.

Each output line has 11 tab-separated fields: sample 1, haplotype 1, sample
2, haplotype 2, chromosome, start and end position (the run's first and last
marker), length in cM (at most 3 decimals), the number of markers, the number
of haplotypes in the VCF file carrying the run's alleles (the pair's own two
included), and that number divided by all haplotypes in the VCF file
(6 significant digits).  As in Refined IBD output, sample 1 is the sample
whose haplotype comes first in the VCF file, and a sample paired with itself
has haplotype 1 first.  A summary goes to standard error.
"""

from __future__ import annotations

import math
import sys
from decimal import ROUND_HALF_EVEN, Decimal, localcontext
from typing import List, NamedTuple, Optional, Set, Tuple

from .genetic_map import PlinkGenMap
from .ids import Indexer
from .input_it import from_gzip_file
from .java_compat import (JavaException, JavaTermination,
                          java_decimal_format, java_encode,
                          java_parse_double, to_java_chars, utilities_exit)
from .pbwt import long_matches
from .vcf import VcfHeader, VcfRecord, read_vcf

PROGRAM = "python -m find_ibs_segments"
DEFAULT_MIN_CM = 3.0


def usage() -> str:
    nl = "\n"
    return (nl
            + "usage: " + PROGRAM + " [vcf] [map] [min_cm] > [out]" + nl
            + nl
            + "where" + nl
            + "  [vcf]    = phased VCF file with biallelic markers (the file "
              "given to Refined IBD)" + nl
            + "  [map]    = PLINK format genetic map file with centimorgan "
              "(cM) distances" + nl
            + "  [min_cm] = optional: minimum run length in cM (default 3; "
              "Refined IBD's length)" + nl
            + "  [out]    = one line per run of identical alleles between "
              "two haplotypes" + nl
            + nl
            + "Every pair of haplotypes, including a sample's own two, is "
              "compared over each" + nl
            + "whole chromosome.  Each run of markers with the same allele "
              "on both haplotypes," + nl
            + "extended until the next marker differs, is written if it is "
              "at least [min_cm]" + nl
            + "cM long, with its number of markers, the number of "
              "haplotypes in [vcf]" + nl
            + "carrying its alleles, and that number divided by all "
              "haplotypes in [vcf]." + nl)


def exit_with_error(message: str) -> JavaTermination:
    """The usage message and the error message on standard error, exit
    status 1."""
    return JavaTermination(stderr=usage() + "\n" + message + "\n", status=1)


def parse_min_cm(arg: str) -> float:
    try:
        min_cm = java_parse_double(arg)
    except JavaException as error:
        if error.java_class == "java.lang.NumberFormatException":
            raise exit_with_error("ERROR: [min_cm] is not a numerical value: "
                                  + arg)
        raise
    if not math.isfinite(min_cm) or min_cm <= 0.0:
        raise exit_with_error("ERROR: [min_cm]=" + arg)
    return min_cm


def format_frequency(carriers: int, n_haplotypes: int) -> str:
    """``carriers / n_haplotypes`` rounded half to even to 6 significant
    digits, in plain decimal notation (``0.0208333``, ``0.5``, ``1``)."""
    with localcontext() as ctx:
        ctx.prec = 6
        ctx.rounding = ROUND_HALF_EVEN
        value = Decimal(carriers) / Decimal(n_haplotypes)
    return format(value, "f")


class Block(NamedTuple):
    """One chromosome's records: ``recs[start:end]``."""
    chrom_index: int
    start: int
    end: int


def chromosome_blocks(recs: List[VcfRecord], chrom_ids: Indexer) -> List[Block]:
    """Each chromosome's records, which must form one block in position
    order (runs are found over consecutive records)."""
    blocks: List[Block] = []
    seen: Set[int] = set()
    for index, rec in enumerate(recs):
        if blocks and rec.chrom_index == blocks[-1].chrom_index:
            prev = recs[index - 1]
            if rec.pos < prev.pos:
                raise utilities_exit(
                    "ERROR: VCF records are not in position order: "
                    "chromosome " + chrom_ids.id(rec.chrom_index)
                    + " position " + str(rec.pos) + " follows position "
                    + str(prev.pos))
            blocks[-1] = blocks[-1]._replace(end=index + 1)
        else:
            if rec.chrom_index in seen:
                raise utilities_exit(
                    "ERROR: the records of chromosome "
                    + chrom_ids.id(rec.chrom_index)
                    + " are not in one block in the VCF file")
            seen.add(rec.chrom_index)
            blocks.append(Block(rec.chrom_index, index, index + 1))
    return blocks


def sample_names(header: VcfHeader, sample_ids: Indexer) -> List[str]:
    """The VCF samples' identifiers, in VCF column order."""
    names = [""] * header.n_samples
    for id_index, column in enumerate(header.id_index_to_index):
        if column != -1:
            names[column] = sample_ids.id(id_index)
    return names


class Result(NamedTuple):
    output: str
    summary: str


def find_runs(recs: List[VcfRecord], block: Block, gen_map: PlinkGenMap,
              chrom_ids: Indexer, n_samples: int,
              min_cm: float) -> List[Tuple[int, int, int, int, int, float]]:
    """The runs on one chromosome: ``(h1, h2, start, end, carriers,
    length)`` with record indices relative to the block."""
    chrom = block.chrom_index
    block_recs = recs[block.start:block.end]
    n_haps = 2 * n_samples
    sites: List[bytes] = []
    for rec in block_recs:
        row = bytearray(n_haps)
        allele1 = rec.allele1
        allele2 = rec.allele2
        for s in range(n_samples):
            a1 = allele1[s]
            a2 = allele2[s]
            if a1 > 1 or a2 > 1:
                raise utilities_exit(
                    "ERROR: find-ibs-segments needs biallelic markers: "
                    "chromosome " + chrom_ids.id(chrom) + " position "
                    + str(rec.pos) + " has allele index " + str(max(a1, a2)))
            row[2 * s] = a1
            row[2 * s + 1] = a2
        sites.append(bytes(row))
    cm = [gen_map.gen_pos(chrom, rec.pos) for rec in block_recs]
    haps = [bytes(site[h] for site in sites) for h in range(n_haps)]
    runs = []
    for h1, h2, start, end in long_matches(sites, cm, min_cm):
        alleles = haps[h1][start:end + 1]
        carriers = sum(1 for hap in haps if hap[start:end + 1] == alleles)
        runs.append((h1, h2, start, end, carriers, cm[end] - cm[start]))
    return runs


def run(args: List[str]) -> Result:
    """Runs the program and returns its output and summary.  Raises
    JavaTermination or JavaException where the program stops."""
    if not args:
        raise JavaTermination(stdout=usage() + "\n", status=0)
    if len(args) not in (2, 3):
        raise exit_with_error("ERROR: expected 2 or 3 arguments, found "
                              + str(len(args)))
    vcf_path, map_path = args[0], args[1]
    min_cm = parse_min_cm(args[2]) if len(args) == 3 else DEFAULT_MIN_CM
    chrom_ids = Indexer()
    sample_ids = Indexer()
    gen_map = PlinkGenMap.from_plink_map_file(map_path, chrom_ids)
    header, recs = read_vcf(vcf_path, chrom_ids, sample_ids,
                            from_gzip_file(vcf_path))
    names = sample_names(header, sample_ids)
    n_haps = 2 * header.n_samples
    blocks = chromosome_blocks(recs, chrom_ids)
    lines: List[Tuple[Tuple, str]] = []
    n_within = 0
    for block_number, block in enumerate(blocks):
        chrom = chrom_ids.id(block.chrom_index)
        for h1, h2, start, end, carriers, length in find_runs(
                recs, block, gen_map, chrom_ids, header.n_samples, min_cm):
            if h1 // 2 == h2 // 2:
                n_within += 1
            line = "\t".join((
                names[h1 // 2], str(h1 % 2 + 1), names[h2 // 2],
                str(h2 % 2 + 1), chrom, str(recs[block.start + start].pos),
                str(recs[block.start + end].pos),
                java_decimal_format(length, 3), str(end - start + 1),
                str(carriers), format_frequency(carriers, n_haps))) + "\n"
            lines.append(((block_number, start, end, h1, h2), line))
    lines.sort(key=lambda item: item[0])
    summary = ("find-ibs-segments: samples " + str(header.n_samples)
               + ", markers " + str(len(recs)) + ", chromosomes "
               + str(len(blocks)) + "; runs of at least "
               + java_decimal_format(min_cm, 3) + " cM: " + str(len(lines))
               + " (" + str(n_within)
               + " between a sample's own two haplotypes)\n")
    return Result("".join(line for _, line in lines), summary)


def main(argv: Optional[List[str]] = None) -> int:
    """Command-line entry point; returns the exit status."""
    args = [to_java_chars(a) for a in (sys.argv[1:] if argv is None
                                       else argv)]
    stdout = sys.stdout.buffer
    stderr = sys.stderr.buffer
    try:
        result = run(args)
        stdout.write(java_encode(result.output))
        stdout.flush()
    except JavaTermination as termination:
        stdout.write(java_encode(termination.stdout))
        stdout.flush()
        stderr.write(java_encode(termination.stderr))
        stderr.flush()
        return termination.status
    except JavaException as error:
        stderr.write(java_encode('Exception in thread "main" '
                                 + error.to_string() + "\n"))
        stderr.flush()
        return 1
    stderr.write(java_encode(result.summary))
    stderr.flush()
    return 0
