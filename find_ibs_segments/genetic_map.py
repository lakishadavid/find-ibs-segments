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
# Python from src/vcf/PlinkGenMap.java in refined-ibd.17Jan20.102.zip (SHA-256
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
"""PLINK genetic map, converted from vcf.PlinkGenMap.

A PLINK map line has four white-space separated fields: chromosome,
identifier, genetic position (cM) and base position.  ``gen_pos`` returns
the genetic position of a base position: the map value at a map position,
linear interpolation between map positions, and, beyond either end of the
map, linear extrapolation along the line through the end position and the
map position about 5 cM inside it.
"""

from __future__ import annotations

import math
from typing import List

from .ids import Indexer
from .input_it import from_gzip_file
from .java_compat import (JavaException, get_fields_ws, get_fields_ws_limit,
                          illegal_argument, java_binary_search_double,
                          java_binary_search_int, java_parse_double,
                          java_parse_int, to_int32)

_MIN_END_CM_DIST = 5.0


def _java_double_to_string_nonfinite(x: float) -> str:
    """``Double.toString`` for a value that is not finite."""
    if x != x:
        return "NaN"
    return "Infinity" if x > 0 else "-Infinity"


class PlinkGenMap:
    """Genetic map positions for each chromosome index."""

    def __init__(self, base_pos: List[List[int]], gen_pos: List[List[float]],
                 chrom_ids: Indexer):
        self._base_pos = base_pos
        self._gen_pos = gen_pos
        self._chrom_ids = chrom_ids

    @classmethod
    def from_plink_map_file(cls, path: str,
                            chrom_ids: Indexer) -> "PlinkGenMap":
        """``PlinkGenMap.fromPlinkMapFile(File)``."""
        chrom_list: List[List[str]] = []
        for line in from_gzip_file(path):
            fields = get_fields_ws_limit(line, 4)
            if len(fields) > 0:
                if len(fields) < 4:
                    raise illegal_argument("Map file format error: " + line)
                chrom_index = chrom_ids.get_index(fields[0])
                while chrom_index >= len(chrom_list):
                    chrom_list.append([])
                chrom_list[chrom_index].append(line)
        base_pos: List[List[int]] = []
        gen_pos: List[List[float]] = []
        for lines in chrom_list:
            b, g = _fill_map_positions(lines)
            base_pos.append(b)
            gen_pos.append(g)
        return cls(base_pos, gen_pos, chrom_ids)

    def _check_chrom_index(self, chrom: int) -> None:
        if chrom < 0 or chrom >= self._chrom_ids.size():
            raise JavaException("java.lang.IndexOutOfBoundsException",
                                str(chrom))
        if chrom >= len(self._base_pos) or len(self._base_pos[chrom]) == 0:
            raise illegal_argument("missing genetic map for chromosome "
                                   + self._chrom_ids.id(chrom))

    def gen_pos(self, chrom: int, base_position: int) -> float:
        """``PlinkGenMap.genPos(int, int)``."""
        self._check_chrom_index(chrom)
        base = self._base_pos[chrom]
        gen = self._gen_pos[chrom]
        map_size_m1 = len(base) - 1
        index = java_binary_search_int(base, base_position)
        if index >= 0:
            return gen[index]
        ins_pt = -index - 1
        a_index = ins_pt - 1
        b_index = ins_pt
        if a_index == map_size_m1:
            ins_pt = java_binary_search_double(
                gen, gen[map_size_m1] - _MIN_END_CM_DIST)
            if ins_pt < 0:
                ins_pt = -ins_pt - 2
            a_index = max(ins_pt, 0)
            b_index = map_size_m1
        elif b_index == 0:
            ins_pt = java_binary_search_double(
                gen, gen[0] + _MIN_END_CM_DIST)
            if ins_pt < 0:
                ins_pt = -ins_pt - 1
            a_index = 0
            b_index = min(ins_pt, map_size_m1)
        x = base_position
        a = base[a_index]
        b = base[b_index]
        fa = gen[a_index]
        fb = gen[b_index]
        return fa + ((float(to_int32(x - a)) / float(to_int32(b - a)))
                     * (fb - fa))


def _fill_map_positions(lines: List[str]):
    """``PlinkGenMap.fillMapPositions``: parses and checks one chromosome's
    map lines."""
    n = len(lines)
    base = [0] * n
    gen = [0.0] * n
    for j, line in enumerate(lines):
        fields = get_fields_ws(line)
        if len(fields) != 4:
            raise illegal_argument("Map file format error: " + line)
        base[j] = java_parse_int(fields[3])
        gen[j] = java_parse_double(fields[2])
        if not math.isfinite(gen[j]):
            raise illegal_argument("invalid map position: "
                                   + _java_double_to_string_nonfinite(gen[j]))
        if j > 0:
            if base[j] == base[j - 1]:
                raise illegal_argument("duplication position: " + line)
            if base[j] < base[j - 1] or gen[j] < gen[j - 1]:
                raise illegal_argument(
                    "map positions not in ascending order: " + line)
    if n > 0 and gen[0] == gen[n - 1]:
        raise illegal_argument("Genetic map has only one map position: "
                               + lines[0])
    return base, gen
