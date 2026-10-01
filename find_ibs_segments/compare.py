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
"""Compares find-ibs-segments runs with Refined IBD's output for the same
VCF file, genetic map and minimum length.

Usage: ``python -m find_ibs_segments.compare [runs] [ibd] [hbd] [log] >
[out]``.

``[runs]`` is find-ibs-segments output; ``[ibd]`` and ``[hbd]`` are Refined
IBD's ``.ibd`` and ``.hbd`` files as Refined IBD wrote them, before any
merging (plain text, or gzip when the name ends in ``.gz``); ``[log]`` is
Refined IBD's ``.log`` file, whose ``Window`` lines give each analysis
window's first and last marker position.

Every Refined IBD segment lies inside one run of the same samples and
haplotypes, because both programs extend a segment over consecutive markers
at which the two haplotypes carry the same allele.  Each run is classified
by the Refined IBD segments inside it:

- ``whole``: one segment with the run's start and end;
- ``start_cut``, ``end_cut``, ``both_cut``: one segment, starting later
  and/or ending earlier than the run;
- ``split``: two or more segments;
- ``not_reported``: none.

A Refined IBD segment inside no run is listed as ``refined_only``; that
would mean the two programs do not see the same markers, or that Refined IBD
does not work as described.  For every segment end inside a run, the output
says whether it is the first or last marker of a Refined IBD window.  A
summary goes to standard error.
"""

from __future__ import annotations

import bisect
import re
import sys
from collections import Counter, defaultdict
from typing import Dict, List, NamedTuple, Optional, Tuple

from .input_it import from_gzip_file
from .java_compat import (JavaException, JavaTermination, java_encode,
                          java_trim, to_java_chars, utilities_exit)

PROGRAM = "python -m find_ibs_segments.compare"

CATEGORIES = ("whole", "start_cut", "end_cut", "both_cut", "split",
              "not_reported", "refined_only")
COLUMNS = ("category", "sample1", "hap1", "sample2", "hap2", "chromosome",
           "run_start", "run_end", "run_cm", "run_markers", "carriers",
           "frequency", "refined_segments", "inner_ends", "windows")

_WINDOW = re.compile(r"Window (\d+) \((?:(.+?):)?(-?\d+)-(?:(.+?):)?(-?\d+)\)")

Key = Tuple[str, str, str, str, str]


def usage() -> str:
    nl = "\n"
    return (nl
            + "usage: " + PROGRAM + " [runs] [ibd] [hbd] [log] > [out]" + nl
            + nl
            + "where" + nl
            + "  [runs] = find-ibs-segments output" + nl
            + "  [ibd]  = Refined IBD .ibd file, before merging (.gz allowed)"
            + nl
            + "  [hbd]  = Refined IBD .hbd file, before merging (.gz allowed)"
            + nl
            + "  [log]  = Refined IBD .log file" + nl
            + "  [out]  = one line per run and per Refined IBD segment inside "
              "no run" + nl)


class Run(NamedTuple):
    key: Key
    start: int
    end: int
    fields: Tuple[str, ...]     # run_cm, run_markers, carriers, frequency


class Window(NamedTuple):
    number: int
    first: int
    last: int


def _lines(path: str):
    for raw in from_gzip_file(path):
        line = java_trim(raw)
        if line:
            yield line


def _int(text: str, path: str, line: str) -> int:
    try:
        return int(text)
    except ValueError:
        raise utilities_exit("ERROR: not an integer position in " + path
                             + ": " + line)


def read_runs(path: str) -> Tuple[Dict[Key, List[Run]], List[str]]:
    """The runs grouped by samples, haplotypes and chromosome, sorted by
    start, and the chromosomes in order of first appearance."""
    runs: Dict[Key, List[Run]] = defaultdict(list)
    chroms: List[str] = []
    for line in _lines(path):
        fields = line.split("\t")
        if len(fields) != 11:
            raise utilities_exit("ERROR: a find-ibs-segments line has 11 "
                                 "fields: " + line)
        key = tuple(fields[:5])
        if fields[4] not in chroms:
            chroms.append(fields[4])
        runs[key].append(Run(key, _int(fields[5], path, line),
                             _int(fields[6], path, line), tuple(fields[7:11])))
    for group in runs.values():
        group.sort(key=lambda r: r.start)
        for a, b in zip(group, group[1:]):
            if b.start <= a.end:
                raise utilities_exit("ERROR: runs overlap in " + path + ": "
                                     + str(a.start) + "-" + str(a.end)
                                     + " and " + str(b.start) + "-"
                                     + str(b.end) + " of "
                                     + " ".join(a.key))
    return runs, chroms


def read_refined(path: str) -> List[Tuple[Key, int, int]]:
    """Refined IBD segments: 9 white-space separated fields."""
    segments = []
    for line in _lines(path):
        fields = line.split()
        if len(fields) != 9:
            raise utilities_exit("ERROR: a Refined IBD line has 9 fields: "
                                 + line)
        segments.append((tuple(fields[:5]), _int(fields[5], path, line),
                         _int(fields[6], path, line)))
    return segments


def read_windows(path: str) -> Dict[str, List[Window]]:
    """Refined IBD's analysis windows, by chromosome, from the ``Window``
    lines of its log.  A window on a chromosome named ``.`` is written
    without a chromosome and is filed under ``.``."""
    windows: Dict[str, List[Window]] = defaultdict(list)
    for line in _lines(path):
        match = _WINDOW.fullmatch(line)
        if match is not None:
            number, chrom, first, _, last = match.groups()
            windows[chrom or "."].append(Window(int(number), int(first),
                                                int(last)))
    if not windows:
        raise utilities_exit("ERROR: no Window lines in " + path)
    return windows


def inner_ends(run: Run, inside: List[Tuple[int, int]],
               windows: List[Window]) -> str:
    """Each segment start or end inside the run (not at the run's own ends),
    with the window edge it falls on, if any."""
    entries = []
    for start, end in inside:
        if start != run.start:
            edge = [str(w.number) for w in windows if w.first == start]
            entries.append("start " + str(start) + " = "
                           + ("first marker of window " + ",".join(edge)
                              if edge else "no window edge"))
        if end != run.end:
            edge = [str(w.number) for w in windows if w.last == end]
            entries.append("end " + str(end) + " = "
                           + ("last marker of window " + ",".join(edge)
                              if edge else "no window edge"))
    return "; ".join(entries) if entries else "."


def classify(run: Run, inside: List[Tuple[int, int]]) -> str:
    if not inside:
        return "not_reported"
    if len(inside) > 1:
        return "split"
    start, end = inside[0]
    if start == run.start and end == run.end:
        return "whole"
    if end == run.end:
        return "start_cut"
    if start == run.start:
        return "end_cut"
    return "both_cut"


class Result(NamedTuple):
    output: str
    summary: str


def run(args: List[str]) -> Result:
    if not args:
        raise JavaTermination(stdout=usage() + "\n", status=0)
    if len(args) != 4:
        raise JavaTermination(stderr=usage() + "\nERROR: expected 4 "
                              "arguments, found " + str(len(args)) + "\n",
                              status=1)
    runs_path, ibd_path, hbd_path, log_path = args
    runs, chroms = read_runs(runs_path)
    windows = read_windows(log_path)
    inside: Dict[Run, List[Tuple[int, int]]] = defaultdict(list)
    refined_only: List[Tuple[Key, int, int]] = []
    starts = {key: [r.start for r in group] for key, group in runs.items()}
    for path in (ibd_path, hbd_path):
        for key, start, end in read_refined(path):
            group = runs.get(key)
            if group is not None:
                i = bisect.bisect_right(starts[key], start) - 1
                if i >= 0 and group[i].end >= end:
                    inside[group[i]].append((start, end))
                    continue
            refined_only.append((key, start, end))
    chrom_order = {c: i for i, c in enumerate(chroms)}
    rows: List[Tuple[Tuple, str]] = []
    counts: Counter = Counter()
    for group in runs.values():
        for r in group:
            segs = sorted(inside.get(r, []))
            category = classify(r, segs)
            counts[category] += 1
            chrom = r.key[4]
            touched = [str(w.number) for w in windows.get(chrom, [])
                       if w.first <= r.end and w.last >= r.start]
            fields = ((category,) + r.key + (str(r.start), str(r.end))
                      + r.fields
                      + (",".join("%d-%d" % s for s in segs) or ".",
                         inner_ends(r, segs, windows.get(chrom, [])),
                         ",".join(touched) or "."))
            rows.append(((chrom_order[chrom], r.start, r.end, r.key),
                         "\t".join(fields) + "\n"))
    for key, start, end in refined_only:
        counts["refined_only"] += 1
        fields = (("refined_only",) + key + (".",) * 6
                  + ("%d-%d" % (start, end), ".", "."))
        rows.append(((chrom_order.get(key[4], len(chrom_order)), start, end,
                      key), "\t".join(fields) + "\n"))
    rows.sort(key=lambda item: item[0])
    output = "\t".join(COLUMNS) + "\n" + "".join(line for _, line in rows)
    summary = ("find-ibs-segments compare: "
               + ", ".join(c + " " + str(counts[c]) for c in CATEGORIES)
               + "\n")
    return Result(output, summary)


def main(argv: Optional[List[str]] = None) -> int:
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


if __name__ == "__main__":
    sys.exit(main())
