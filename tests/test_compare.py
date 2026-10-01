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
"""Tests of the comparison with Refined IBD output.

Run from the repository root: ``python -m unittest discover -s tests``.

The runs come from the finder on the hand-built data of test_find.py; the
Refined IBD files are written by hand to give each category once.  All data
are synthetic.
"""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from find_ibs_segments import compare, find  # noqa: E402
from find_ibs_segments.java_compat import JavaTermination  # noqa: E402
from test_find import GENO, write_inputs  # noqa: E402

IBD = ["A\t1\tB\t1\t1\t1000\t20000\t5.1\t2.375",       # whole
       "A\t1\tB\t2\t1\t5000\t20000\t4.2\t1.875",       # start cut
       "A\t2\tB\t2\t1\t1000\t30000\t6.0\t3.625",       # end cut at window 1
       "A\t2\tB\t1\t1\t1000\t15000\t2.0\t1.75",        # split ...
       "A\t2\tB\t1\t1\t25000\t40000\t2.0\t1.875",      # ... at window 2
       "A\t1\tB\t2\t1\t23000\t33000\t1.0\t1.25",       # cut at both ends
       "A\t1\tB\t1\t1\t22000\t40000\t3.3\t2.25",       # whole
       "A\t1\tB\t1\t1\t18000\t25000\t0.5\t0.875"]      # across a mismatch
HBD = ["A\t1\tA\t2\t1\t1000\t20000\t7.0\t2.375"]       # whole
LOG = ["", "Window 1 (1:1000-30000)", "Markers: 30", "",
       "Window 2 (1:25000-40000)", "Markers: 16"]
HEADER = "\t".join(compare.COLUMNS)


def line(category, key, start, end, run_fields, refined, inner, windows):
    return "\t".join((category,) + key + (str(start), str(end)) + run_fields
                     + (refined, inner, windows))


EXPECTED = [
    line("whole", ("A", "1", "A", "2", "1"), 1000, 20000, ("2.375", "20", "4", "1"),
         "1000-20000", ".", "1"),
    line("whole", ("A", "1", "B", "1", "1"), 1000, 20000, ("2.375", "20", "4", "1"),
         "1000-20000", ".", "1"),
    line("start_cut", ("A", "1", "B", "2", "1"), 1000, 20000,
         ("2.375", "20", "4", "1"), "5000-20000",
         "start 5000 = no window edge", "1"),
    line("end_cut", ("A", "2", "B", "2", "1"), 1000, 34000,
         ("4.125", "34", "3", "0.75"), "1000-30000",
         "end 30000 = last marker of window 1", "1,2"),
    line("not_reported", ("B", "1", "B", "2", "1"), 1000, 34000,
         ("4.125", "34", "3", "0.75"), ".", ".", "1,2"),
    line("split", ("A", "2", "B", "1", "1"), 1000, 40000,
         ("4.875", "40", "2", "0.5"), "1000-15000,25000-40000",
         "end 15000 = no window edge; start 25000 = first marker of window 2",
         "1,2"),
    "\t".join(("refined_only", "A", "1", "B", "1", "1") + (".",) * 6
              + ("18000-25000", ".", ".")),
    line("both_cut", ("A", "1", "B", "2", "1"), 22000, 34000,
         ("1.5", "13", "4", "1"), "23000-33000",
         "start 23000 = no window edge; end 33000 = no window edge", "1,2"),
    line("not_reported", ("A", "1", "A", "2", "1"), 22000, 40000,
         ("2.25", "19", "3", "0.75"), ".", ".", "1,2"),
    line("whole", ("A", "1", "B", "1", "1"), 22000, 40000,
         ("2.25", "19", "3", "0.75"), "22000-40000", ".", "1,2"),
]


def files(d, ibd=IBD, hbd=HBD, log=LOG, gz=True):
    vcf, mp = write_inputs(d, GENO)
    runs = os.path.join(d, "runs.ibs")
    with open(runs, "wb") as f:
        f.write(find.run([vcf, mp, "1.5"]).output.encode())
    paths = []
    for name, rows in (("out.ibd", ibd), ("out.hbd", hbd)):
        path = os.path.join(d, name + (".gz" if gz else ""))
        opener = gzip.open if gz else open
        with opener(path, "wb") as f:
            f.write("".join(r + "\n" for r in rows).encode())
        paths.append(path)
    log_path = os.path.join(d, "out.log")
    with open(log_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(log) + "\n")
    return [runs] + paths + [log_path]


class Compare(unittest.TestCase):

    def test_every_category(self):
        for gz in (True, False):
            with self.subTest(gz=gz):
                with tempfile.TemporaryDirectory() as d:
                    got = compare.run(files(d, gz=gz))
                self.assertEqual(got.output.splitlines(), [HEADER] + EXPECTED)
                self.assertEqual(got.summary,
                                 "find-ibs-segments compare: whole 3, start_cut 1, "
                                 "end_cut 1, both_cut 1, split 1, not_reported 2, "
                                 "refined_only 1\n")

    def test_command_line(self):
        with tempfile.TemporaryDirectory() as d:
            args = files(d)
            done = subprocess.run([sys.executable, "-m", "find_ibs_segments.compare"]
                                  + args, capture_output=True, cwd=REPO)
        self.assertEqual(done.returncode, 0)
        self.assertEqual(done.stdout.decode().splitlines(), [HEADER] + EXPECTED)

    def test_errors(self):
        cases = (
            (dict(ibd=["A\t1\tB\t1\t1\t1000\t20000\t5.1"]),
             "ERROR: a Refined IBD line has 9 fields: A\t1\tB\t1\t1\t1000\t20000\t5.1\n"),
            (dict(log=["Markers: 30"]), "no Window lines"),
        )
        for kw, message in cases:
            with self.subTest(message=message):
                with tempfile.TemporaryDirectory() as d:
                    with self.assertRaises(JavaTermination) as ctx:
                        compare.run(files(d, **kw))
                self.assertIn(message, ctx.exception.stderr)
        with self.assertRaises(JavaTermination) as ctx:
            compare.run(["a", "b"])
        self.assertEqual(ctx.exception.status, 1)


if __name__ == "__main__":
    unittest.main()
