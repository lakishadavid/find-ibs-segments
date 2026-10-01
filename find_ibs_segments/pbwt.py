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
"""Maximal identical runs between haplotypes, found with the positional
Burrows-Wheeler transform (PBWT).

Reference: Durbin R (2014).  Efficient haplotype matching and storage using
the positional Burrows-Wheeler transform (PBWT).  Bioinformatics
30(9):1266-1272.  doi:10.1093/bioinformatics/btu014.  The prefix and
divergence arrays are updated as in the paper's Algorithm 2; long matches are
reported as in its Algorithm 3, with the length measured in cM instead of
markers.  Written from the paper; no code was copied.
"""

from __future__ import annotations

from typing import Iterator, List, Sequence, Tuple


def long_matches(sites: Sequence[bytes], cm: Sequence[float],
                 min_cm: float) -> Iterator[Tuple[int, int, int, int]]:
    """Every maximal run of markers over which two haplotypes carry the same
    allele, with ``cm[end] - cm[start] >= min_cm``.

    ``sites[k][h]`` is haplotype ``h``'s allele (0 or 1) at marker ``k``;
    ``cm`` holds the markers' genetic positions, in non-decreasing order;
    ``min_cm`` must be positive.  Yields ``(h1, h2, start, end)`` with
    ``h1 < h2`` and ``start``/``end`` the first and last marker of the run,
    so that marker ``start - 1`` (if any) and marker ``end + 1`` (if any)
    carry different alleles on the two haplotypes.  Each run is yielded
    once; the order is not specified.
    """
    n_sites = len(sites)
    if n_sites == 0:
        return
    n_haps = len(sites[0])
    # a: haplotypes sorted by their reversed prefixes up to the current
    # marker; d[i]: the first marker of the match between a[i] and a[i - 1]
    # ending just before the current marker (Durbin's divergence array).
    a: List[int] = list(range(n_haps))
    d: List[int] = [0] * n_haps
    first_long = 0
    for k in range(n_sites + 1):
        if k > 0:
            # Matches ending at marker k - 1 are long when they start at or
            # before marker m.
            while first_long < k and cm[k - 1] - cm[first_long] >= min_cm:
                first_long += 1
            m = first_long - 1
            alleles = sites[k] if k < n_sites else None
            i = 0
            while i < n_haps:
                j = i + 1
                while j < n_haps and d[j] <= m:
                    j += 1
                # a[i:j] is a block: every pair in it matches over at least
                # [m, k - 1].  A pair whose alleles differ at marker k (or
                # every pair, after the last marker) ends its run at k - 1.
                if j - i > 1 and (alleles is None or _mixed(alleles, a, i, j)):
                    for q in range(i + 1, j):
                        hq = a[q]
                        start = 0
                        for p in range(q - 1, i - 1, -1):
                            if d[p + 1] > start:
                                start = d[p + 1]
                            hp = a[p]
                            if alleles is None or alleles[hp] != alleles[hq]:
                                if hp < hq:
                                    yield hp, hq, start, k - 1
                                else:
                                    yield hq, hp, start, k - 1
                i = j
        if k == n_sites:
            break
        # Durbin's Algorithm 2: sort by allele at marker k, keeping order.
        site = sites[k]
        p = q = k + 1
        a0: List[int] = []
        a1: List[int] = []
        d0: List[int] = []
        d1: List[int] = []
        for i in range(n_haps):
            h = a[i]
            di = d[i]
            if di > p:
                p = di
            if di > q:
                q = di
            if site[h] == 0:
                a0.append(h)
                d0.append(p)
                p = 0
            else:
                a1.append(h)
                d1.append(q)
                q = 0
        a = a0 + a1
        d = d0 + d1


def _mixed(alleles: bytes, a: List[int], i: int, j: int) -> bool:
    """Whether the haplotypes ``a[i:j]`` carry both alleles."""
    first = alleles[a[i]]
    return any(alleles[a[x]] != first for x in range(i + 1, j))
