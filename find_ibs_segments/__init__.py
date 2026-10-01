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
"""find-ibs-segments: runs of identical alleles between haplotypes of a
phased VCF file, the first part of an IBS segment detector.

Modules: ``find`` (the program), ``pbwt`` (positional Burrows-Wheeler
transform), ``compare`` (comparison with Refined IBD output), and, copied
from merge-ibs-segments, ``vcf`` (VCF reading), ``genetic_map`` (PLINK map),
``ids`` (identifier indexes), ``input_it`` (line input) and ``java_compat``
(Java behaviour those depend on).
"""
