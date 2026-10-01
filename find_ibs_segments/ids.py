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
# Python from src/beagleutil/ChromIds.java, src/beagleutil/SampleIds.java and
# src/beagleutil/ThreadSafeIndexer.java in refined-ibd.17Jan20.102.zip (SHA-256
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
"""Identifier indexes, converted from beagleutil.ChromIds, SampleIds and
ThreadSafeIndexer.

In the Java program ChromIds and SampleIds are global: an identifier gets the
next index the first time any part of the program asks for it.  The genetic
map registers chromosomes first, the VCF file registers its samples (in
column order) and its chromosomes next, and the IBD segments register the
remaining identifiers.  The Java program prints chromosome indices, not names,
in its "position missing" message, so the order of registration is kept.
"""

from __future__ import annotations

from typing import Dict, List

from .java_compat import illegal_argument


class Indexer:
    """An index of string identifiers: each distinct identifier gets the next
    index, starting at 0."""

    def __init__(self) -> None:
        self._index: Dict[str, int] = {}
        self._ids: List[str] = []

    def get_index(self, identifier: str) -> int:
        """``getIndex(String)``: the identifier's index, registering it if
        it is new."""
        if not identifier:
            raise illegal_argument("id.isEmpty()")
        index = self._index.get(identifier)
        if index is None:
            index = len(self._ids)
            self._ids.append(identifier)
            self._index[identifier] = index
        return index

    def get_indices(self, identifiers: List[str]) -> List[int]:
        """``getIndices(String[])``: checks every identifier for emptiness
        first, then registers them in order."""
        for identifier in identifiers:
            if not identifier:
                raise illegal_argument("id.isEmpty()")
        return [self.get_index(identifier) for identifier in identifiers]

    def id(self, index: int) -> str:
        """``id(int)``: the identifier with this index."""
        return self._ids[index]

    def size(self) -> int:
        """``size()``: the number of identifiers."""
        return len(self._ids)
