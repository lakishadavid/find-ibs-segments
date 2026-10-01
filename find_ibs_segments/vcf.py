# find-ibs-segments.  This module comes from merge-ibs-segments, a modified
# Python translation of merge-ibd-segments from Refined-IBD version 17Jan20.102.
#
# Copyright (C) 2014-2016 Brian L. Browning
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
#
# Modified by LaKisha T. David, September 29-30, 2026: translated from Java to
# Python from src/vcf/VcfHeader.java, src/vcf/VcfMetaInfo.java,
# src/vcf/RefIt.java, src/vcf/VcfRecGTParser.java, src/vcf/BasicMarker.java,
# src/vcf/RefGTRec.java, src/vcf/LowMafRefGT.java,
# src/vcf/LowMafRefDiallelicGT.java, src/beagleutil/Samples.java,
# src/bref/SeqCoder3.java and src/blbutil/Const.java in
# refined-ibd.17Jan20.102.zip (SHA-256
# 0b0abf48528ec53d3fec7f21a9742c6e5960857c5a225b0e49346179bc45e6ea). Later
# changes are recorded in the Git history.
#
# Copied into find-ibs-segments on October 1, 2026, with only the first lines
# of this notice and this paragraph changed, from
# merge-ibs-segments (https://github.com/lakishadavid/merge-ibs-segments),
# commit 93402c9.  Later changes are recorded in the Git history.
#
# Brian L. Browning's original notice, reproduced unchanged from those files:
#
# Copyright (C) 2014-2016 Brian L. Browning
#
# This file is part of Refined-IBD
#
# Refined-IBD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Beagle is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Phased VCF reading, converted from the classes listed in the notice above
(chiefly vcf.VcfHeader, vcf.RefIt, vcf.VcfRecGTParser and vcf.BasicMarker).

merge-ibd-segments reads the VCF file through RefIt, which accepts only
phased, non-missing genotypes.  For each record it keeps the marker
(chromosome and position) and the two alleles of every sample.  The checks
and messages below follow the Java classes line by line.

RefIt reads lines on one thread and parses records on others, so when a file
contains more than one error, which error the Java program reports can
depend on timing.  This conversion reads and parses in file order and reports
the first error in that order.
"""

from __future__ import annotations

import array
import re
from typing import Iterator, List, Optional, Tuple

from .ids import Indexer
from .java_compat import (JavaException, JavaTermination, get_fields_delim,
                          get_fields_delim_limit, illegal_argument,
                          java_array_to_string, java_file_string, java_is_acgtn,
                          java_is_digit, java_is_whitespace, java_parse_int,
                          java_substring, java_trim, utilities_exit,
                          utilities_exit_throwable)

_SHORT_HEADER_PREFIX = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
_HEADER_PREFIX = _SHORT_HEADER_PREFIX + "\tFORMAT"
_SAMPLE_OFFSET = 9
_PURE_GT = re.compile("(?:[0-9]\\|[0-9]\t)*[0-9]\\|[0-9]")
_DIGIT_VALUES = bytes(range(256)).translate(
    bytes.maketrans(b"0123456789", bytes(range(10))))


class VcfRecord:
    """A marker and the phased alleles of every sample.

    ``allele1[s]`` and ``allele2[s]`` are the first and second alleles of
    sample ``s`` (``GTRec.allele1(s)`` and ``allele2(s)``).
    """

    __slots__ = ("chrom_index", "pos", "allele1", "allele2")

    def __init__(self, chrom_index: int, pos: int, allele1, allele2):
        self.chrom_index = chrom_index
        self.pos = pos
        self.allele1 = allele1
        self.allele2 = allele2


class _Marker:
    """The fields of vcf.BasicMarker that the program uses."""

    __slots__ = ("chrom_index", "pos", "ids", "alleles", "end")

    def __init__(self, chrom_index: int, pos: int, ids: List[str],
                 alleles: List[str], end: int):
        self.chrom_index = chrom_index
        self.pos = pos
        self.ids = ids
        self.alleles = alleles
        self.end = end

    def to_string(self, chrom_ids: Indexer) -> str:
        """``BasicMarker.toString()``."""
        parts = [chrom_ids.id(self.chrom_index), "\t", str(self.pos)]
        if not self.ids:
            parts.append("\t.")
        else:
            parts.append("\t" + ";".join(self.ids))
        if len(self.alleles) == 1:
            parts.append("\t" + self.alleles[0] + "\t.")
        else:
            parts.append("\t" + self.alleles[0] + "\t"
                         + ",".join(self.alleles[1:]))
        return "".join(parts)


class VcfHeader:
    """``vcf.VcfHeader``: meta-information lines, header line and samples."""

    def __init__(self, lines: Iterator[str], file_display: str,
                 sample_ids: Indexer):
        candidate: Optional[str] = None
        for raw in lines:
            line = java_trim(raw)
            if line.startswith("##"):
                _check_meta_info(line)
            else:
                candidate = line
                break
        _check_header_line(candidate, file_display)
        header_fields = get_fields_delim(candidate, "\t")
        self.file_display = file_display
        self.n_header_fields = len(header_fields)
        ids = header_fields[_SAMPLE_OFFSET:]
        id_indices = sample_ids.get_indices(ids)
        # beagleutil.Samples: the index of a sample identifier in the VCF.
        id_index_to_index = [-1] * sample_ids.size()
        for j, id_index in enumerate(id_indices):
            if id_index_to_index[id_index] != -1:
                raise illegal_argument(
                    "duplicate sample: " + sample_ids.id(id_index)
                    + " (ID index: " + str(id_index) + ")")
            id_index_to_index[id_index] = j
        self.id_index_to_index = id_index_to_index
        self.n_samples = len(ids)

    def sample_index(self, id_index: int) -> int:
        """``Samples.index(int)``: the VCF column index of a sample
        identifier index, or -1 if the sample is not in the VCF file."""
        if id_index >= len(self.id_index_to_index):
            return -1
        return self.id_index_to_index[id_index]


def _check_meta_info(line: str) -> None:
    """``VcfMetaInfo(String)``."""
    index = line.find("=")
    if index <= 0 or index == len(line) - 1:
        raise illegal_argument('VCF meta-information line: missing "="')


def _check_header_line(line: Optional[str], file_display: str) -> None:
    """``VcfHeader.checkHeaderLine``."""
    if line is None or not line.startswith("#"):
        raise illegal_argument(
            "Missing line (#CHROM ...) after meta-information lines\n"
            "File source: " + file_display + "\n"
            + ("null" if line is None else line))
    if not line.startswith(_HEADER_PREFIX):
        if line != _SHORT_HEADER_PREFIX:
            raise illegal_argument(
                "Missing header line (file source: " + file_display + ")\n"
                "The first line after the initial meta-information lines\n"
                "does not begin with: \n" + _HEADER_PREFIX + "\n" + line
                + "\nThe data fields in the header line must be "
                "tab-separated.")


def _data_lines(lines: Iterator[str]) -> Iterator[str]:
    """``RefIt.readLine`` repeated: skips lines that are empty after
    trimming, except the last line of the file, which is returned even if
    blank."""
    blank: Optional[str] = None
    for line in lines:
        if java_trim(line) == "":
            blank = line
            continue
        blank = None
        yield line
    if blank is not None:
        yield blank


def read_vcf(path: str, chrom_ids: Indexer, sample_ids: Indexer,
             lines: Iterator[str]) -> Tuple[VcfHeader, List[VcfRecord]]:
    """Reads the header and all records, as ``RefIt.create`` and the
    ``MergeIbdSegments.Data`` constructor do."""
    file_display = java_file_string(path)
    header = VcfHeader(lines, file_display, sample_ids)
    if header.n_samples < 1:
        # RefIt's constructor builds a bref.SeqCoder3, whose
        # defaultMaxNSeq(nSamples) rejects a sample count below 1 before
        # any record is read.
        raise illegal_argument(str(header.n_samples))
    records: List[VcfRecord] = []
    first = True
    for line in _data_lines(lines):
        if first:
            first = False
            if len(line) == 0:
                # RefIt.stringBufferSize divides by twice the length of the
                # first line; the reading thread reports the failure through
                # Utilities.exit(Throwable).
                raise utilities_exit_throwable(JavaException(
                    "java.lang.ArithmeticException", "/ by zero"))
        if "\t" not in line:
            raise utilities_exit(
                "\nERROR: Missing tab delimiter in VCV Record:\n" + line
                + "\nExiting Program")
        records.append(_parse_record(line, header, chrom_ids))
    if not records:
        raise illegal_argument(
            "No VCF records found (data source: " + file_display + ")\n"
            "Check that the chromosome identifiers are the same in each "
            "input VCF\n"
            "file and in the 'chrom=' command line argument (if 'chrom=' "
            "is used).")
    return header, records


# ---------------------------------------------------------------------------
# vcf.BasicMarker(String)
# ---------------------------------------------------------------------------

def _coordinate(chrom_ids: Indexer, chrom: int, pos: int) -> str:
    return chrom_ids.id(chrom) + ":" + str(pos)


def _parse_marker(rec: str, chrom_ids: Indexer) -> _Marker:
    fields = get_fields_delim_limit(rec, "\t", 9)
    if len(fields) < 8:
        raise utilities_exit(
            "VCF record does not contain at least 8 tab-delimited fields: "
            + rec)
    chrom_index = _extract_chrom(fields[0], rec, chrom_ids)
    pos = _extract_pos(fields[1], rec)
    ids = _extract_ids(chrom_ids, chrom_index, pos, fields[2])
    alleles = _extract_alleles(chrom_ids, chrom_index, pos, fields[3],
                               fields[4])
    end = _extract_end(chrom_ids, chrom_index, pos, rec)
    return _Marker(chrom_index, pos, ids, alleles, end)


def _extract_chrom(chrom: str, rec: str, chrom_ids: Indexer) -> int:
    if chrom == "" or chrom == ".":
        raise utilities_exit("ERROR: missing CHROM field: \n"
                             + java_substring(rec, 0, 80))
    for c in chrom:
        if c == ":" or java_is_whitespace(c):
            raise utilities_exit("invalid character in CHROM field ['" + c
                                 + "']: \n" + java_substring(rec, 0, 80))
    return chrom_ids.get_index(chrom)


def _extract_pos(pos: str, rec: str) -> int:
    for c in pos:
        if not java_is_digit(c):
            raise utilities_exit("ERROR: invalid POS field [" + pos + "]: \n"
                                 + java_substring(rec, 0, 80))
    return java_parse_int(pos)


def _remove_missing_ids(chrom_ids: Indexer, chrom: int, pos: int,
                        ids: List[str]) -> List[str]:
    """``BasicMarker.removeMissingIds``.  It moves the identifiers it keeps
    to the front of ``ids`` in place; ``extractIds`` ignores its result and
    keeps the modified array."""
    index = 0
    for j in range(len(ids)):
        identifier = ids[j]
        for c in identifier:
            if java_is_whitespace(c):
                raise utilities_exit(
                    "ERROR: ID field contains white-space at "
                    + _coordinate(chrom_ids, chrom, pos) + " [" + identifier
                    + "]")
        if len(identifier) > 0 and identifier != ".":
            ids[index] = identifier
            index += 1
    return ids[:index] if index < len(ids) else ids


def _extract_ids(chrom_ids: Indexer, chrom: int, pos: int,
                 field: str) -> List[str]:
    if field == "":
        raise utilities_exit("ERROR: missing ID field at "
                             + _coordinate(chrom_ids, chrom, pos))
    if field == ".":
        return []
    ids = get_fields_delim(field, ";")
    _remove_missing_ids(chrom_ids, chrom, pos, ids)
    return ids


def _extract_alleles(chrom_ids: Indexer, chrom: int, pos: int, ref: str,
                     alt: str) -> List[str]:
    if ref == "":
        raise utilities_exit("ERROR: missing REF field at "
                             + _coordinate(chrom_ids, chrom, pos))
    if alt == "":
        raise utilities_exit("ERROR: missing ALT field: at "
                             + _coordinate(chrom_ids, chrom, pos))
    alt_alleles = [] if alt == "." else get_fields_delim(alt, ",")
    alleles = [ref] + alt_alleles
    _check_alleles(chrom_ids, chrom, pos, alleles)
    return alleles


def _check_alleles(chrom_ids: Indexer, chrom: int, pos: int,
                   alleles: List[str]) -> None:
    coordinate = _coordinate(chrom_ids, chrom, pos)
    if len(set(alleles)) != len(alleles):
        raise utilities_exit("ERROR: duplicate allele at " + coordinate + " "
                             + java_array_to_string(alleles))
    ref = alleles[0]
    if ref == "":
        raise utilities_exit("ERROR: missing REF field at " + coordinate)
    for c in ref:
        if not java_is_acgtn(c):
            raise utilities_exit(
                "ERROR: REF field is not a sequence of A, C, T, G, or N "
                "characters at " + coordinate + " [" + ref + "]")
    for alt in alleles[1:]:
        _check_alt_allele(coordinate, alt)


def _check_alt_allele(coordinate: str, alt: str) -> None:
    n = len(alt)
    if n == 1 and alt[0] == "*":
        return
    if n >= 2 and alt[0] == "<" and alt[n - 1] == ">":
        for c in alt[1:n - 1]:
            if java_is_whitespace(c) or c in ",<>":
                raise utilities_exit("ERROR: invalid ALT allele at "
                                     + coordinate + " [" + alt + "]")
    else:
        for c in alt:
            if not java_is_acgtn(c):
                raise utilities_exit("ERROR: invalid ALT allele at "
                                     + coordinate + " [" + alt + "]")


def _extract_end(chrom_ids: Indexer, chrom: int, pos: int, rec: str) -> int:
    """``BasicMarker.extractEnd``: Java splits the whole record, not only
    the INFO field, at ``;`` and reads every piece that begins ``END=``."""
    if not rec.startswith("END=") and ";END=" not in rec:
        return -1
    end = -1
    for field in get_fields_delim(rec, ";"):
        if field.startswith("END="):
            value = field[4:]
            for c in value:
                if not java_is_digit(c):
                    raise utilities_exit(
                        "ERROR: invalid INFO:END field at "
                        + _coordinate(chrom_ids, chrom, pos) + " [END="
                        + value + "]")
            end = java_parse_int(value)
            if end != -1 and end < pos:
                raise utilities_exit(
                    "ERROR: invalid INFO:END field at "
                    + _coordinate(chrom_ids, chrom, pos) + " [" + str(end)
                    + "]")
    return end


# ---------------------------------------------------------------------------
# vcf.VcfRecGTParser and RefGTRec.alleleCodedInstance
# ---------------------------------------------------------------------------

def _parse_record(rec: str, header: VcfHeader,
                  chrom_ids: Indexer) -> VcfRecord:
    n_samples = header.n_samples
    marker = _parse_marker(rec, chrom_ids)
    n_alleles = len(marker.alleles)
    ninth = -1
    for _ in range(9):
        ninth = rec.find("\t", ninth + 1)
        if ninth == -1:
            raise illegal_argument("VCF record format error: " + rec)
    alleles = _fast_phased_alleles(rec, ninth, n_samples, n_alleles)
    if alleles is None:
        alleles = _phased_alleles(rec, ninth, n_samples, n_alleles, marker,
                                  header, chrom_ids)
    return VcfRecord(marker.chrom_index, marker.pos, alleles[0], alleles[1])


def _fast_phased_alleles(rec: str, ninth: int, n_samples: int,
                         n_alleles: int):
    """The alleles when every sample's genotype is ``d|d`` with single
    digits (optionally followed by ``:`` and other FORMAT fields), all
    below the number of alleles; otherwise None, and the exact parser runs.
    For such records the exact parser returns the same alleles."""
    region = rec[ninth + 1:]
    need = 4 * n_samples - 1
    if (len(region) >= need
            and (len(region) == need or region[need] == "\t")
            and _PURE_GT.fullmatch(region, 0, need) is not None):
        gts = region[:need]
    else:
        fields = region.split("\t", n_samples)
        if len(fields) < n_samples:
            return None
        gts = "\t".join(f.partition(":")[0] for f in fields[:n_samples])
        if len(gts) != need or _PURE_GT.fullmatch(gts) is None:
            return None
    a1 = gts[0::4].encode("ascii").translate(_DIGIT_VALUES)
    a2 = gts[2::4].encode("ascii").translate(_DIGIT_VALUES)
    if max(a1) >= n_alleles or max(a2) >= n_alleles:
        return None
    return a1, a2


def _phased_alleles(rec: str, ninth: int, n_samples: int, n_alleles: int,
                    marker: _Marker, header: VcfHeader, chrom_ids: Indexer):
    """``VcfRecGTParser.phasedAlleles()``."""
    a1s: List[int] = []
    a2s: List[int] = []
    pos = ninth
    for _ in range(n_samples):
        if pos == -1:
            raise _field_count_error(rec, header)
        al_end1 = _al_end1(rec, pos + 1)
        al_end2 = _al_end2(rec, al_end1 + 1)
        a1 = _parse_allele(rec, pos + 1, al_end1, n_alleles)
        a2 = _parse_allele(rec, al_end1 + 1, al_end2, n_alleles)
        if rec[al_end1] != "|" or a1 == -1 or a2 == -1:
            raise illegal_argument(
                "Unphased or missing reference genotype at marker: "
                + marker.to_string(chrom_ids))
        a1s.append(a1)
        a2s.append(a2)
        pos = rec.find("\t", al_end2)
    if max(max(a1s), max(a2s)) < 256:
        return bytes(a1s), bytes(a2s)
    return array.array("i", a1s), array.array("i", a2s)


def _al_end1(rec: str, start: int) -> int:
    n = len(rec)
    if start == n:
        raise _gt_format_error(rec, n)
    index = start
    while index < n:
        c = rec[index]
        if c == "/" or c == "|":
            return index
        if c == ":" or c == "\t":
            raise _gt_format_error(rec, index + 1)
        index += 1
    raise _gt_format_error(rec, n)


def _al_end2(rec: str, start: int) -> int:
    n = len(rec)
    index = start
    while index < n:
        c = rec[index]
        if c == ":" or c == "\t":
            return index
        index += 1
    return index


def _parse_allele(rec: str, start: int, end: int, n_alleles: int) -> int:
    if start == end:
        raise illegal_argument("Missing sample allele: " + rec)
    if start + 1 == end:
        c = rec[start]
        if c == ".":
            return -1
        allele = ord(c) - 48
    else:
        allele = java_parse_int(rec[start:end])
    if allele < 0 or allele >= n_alleles:
        raise illegal_argument("invalid allele [" + rec[start:end] + "]: \n"
                               + rec)
    return allele


def _gt_format_error(rec: str, index: int) -> JavaTermination:
    """``VcfRecGTParser.throwGTFormatError``."""
    return utilities_exit("ERROR: genotype is missing allele separator:\n"
                          + rec[:index] + "\nExiting Program\n")


def _field_count_error(rec: str, header: VcfHeader) -> JavaTermination:
    """``VcfRecGTParser.throwFieldCountError``."""
    fields = get_fields_delim(rec, "\t")
    return utilities_exit(
        "VCF header line has " + str(header.n_header_fields)
        + " fields, but data line has " + str(len(fields)) + " fields\n"
        "File source: " + header.file_display + "\n"
        + java_array_to_string(fields) + "\n")
