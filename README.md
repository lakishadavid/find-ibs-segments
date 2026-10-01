# find-ibs-segments

Lists every run of identical alleles, at least a given length in cM, between
two haplotypes of a phased VCF file, and compares the runs with Brian L.
Browning's [Refined IBD](https://faculty.washington.edu/browning/refined-ibd.html)
output for the same file.

This is the first part of an IBS segment detector. It finds the runs the
way Refined IBD extends its segments, without Refined IBD's analysis windows
and without its LOD tests. Comparing the two shows which runs Refined IBD
reports whole, cuts short, splits or does not report.

## What it does

Every pair of haplotypes is compared over each whole chromosome. That covers
haplotypes of two samples and also a sample's own two haplotypes. A run is a
stretch of consecutive markers at which the two haplotypes carry the same
allele, extended in both directions until the next marker differs or the
chromosome's markers end. Refined-IBD 17Jan20.102 extends its segments the
same way (`IbsHapSegments.inclStart` and `inclEnd`). The run's length is the
genetic map distance between its first and last marker. A run is written if
that length is at least `[min_cm]`, the same test Refined IBD applies with its
`length` parameter (`IbsHapSegments.java:476`).

The runs are found with the positional Burrows–Wheeler transform (Durbin R,
2014, Efficient haplotype matching and storage using the positional
Burrows-Wheeler transform (PBWT), *Bioinformatics* 30(9):1266-1272,
doi:10.1093/bioinformatics/btu014). Its time grows with the number of
haplotypes times the number of markers, plus the number of runs written,
instead of comparing every pair. The code was written from the paper.

Refined IBD differs in three ways (read from its source):
- **Seeds:** it starts a segment only from a 1.5 cM block that is identical
  throughout, on a block grid that starts at each window's first marker.
- **Windows:** it searches 40 cM windows that overlap by at least 7.5 cM, and
  joins the pieces of a run at a splice marker.
- **LOD tests:** it keeps a piece only if its frequency LOD and its IBD LOD
  pass, computed with each window's own haplotype model.

So a run found here can be missing from Refined IBD's output, cut short at a
window edge, or split. The comparison shows where.

## Use

No installation and no packages beyond the Python standard library are
needed; the code has been run with Python 3.12 and 3.13. From the repository
root:

```
python -m find_ibs_segments phased.vcf.gz plink.map 3 > runs.ibs
python -m find_ibs_segments.compare runs.ibs out.ibd.gz out.hbd.gz out.log > comparison.tsv
```

### Finding runs

| Argument | Meaning |
|---|---|
| `[vcf]` | phased VCF file with biallelic markers (plain text, or gzip/BGZF when the name ends in `.gz`); for a comparison, the file given to Refined IBD |
| `[map]` | PLINK genetic map with cM positions (plain text, or gzip/BGZF) |
| `[min_cm]` | optional: minimum run length in cM (default 3); for a comparison, Refined IBD's `length` |

Each chromosome's records must form one block in position order. A missing
or unphased genotype, or a marker with more than two alleles, stops the
program with an error.

Output: one line per run, 11 tab-separated fields, sorted by chromosome,
start, end and haplotypes:

1. sample 1;
2. haplotype of sample 1 (1 or 2);
3. sample 2;
4. haplotype of sample 2;
5. chromosome;
6. start: position of the run's first marker;
7. end: position of the run's last marker;
8. length in cM, with at most 3 decimals;
9. markers: the number of markers in the run;
10. carriers: the number of haplotypes in the VCF file that carry the run's
    alleles, including the pair's own two;
11. frequency: carriers divided by all haplotypes in the VCF file, rounded
    half to even to 6 significant digits.

As in Refined IBD output, sample 1 is the sample whose haplotype comes first
in the VCF file, and a sample paired with itself has haplotype 1 first. The
carriers and the frequency depend on who is in the VCF file. A summary goes
to standard error.

### Comparing with Refined IBD

| Argument | Meaning |
|---|---|
| `[runs]` | output of the finder |
| `[ibd]`, `[hbd]` | Refined IBD's `.ibd` and `.hbd` files as it wrote them, before any merging (plain text, or gzip when the name ends in `.gz`) |
| `[log]` | Refined IBD's `.log` file; its `Window` lines give each analysis window's first and last marker position |

Run Refined IBD and the finder on the same VCF file and map, with the same
minimum length. Every Refined IBD segment then lies inside one run of the
same samples and haplotypes. The output is tab-separated, with a header line
and one line per run:

| Category | Refined IBD segments inside the run |
|---|---|
| `whole` | one, with the run's start and end |
| `start_cut` / `end_cut` / `both_cut` | one, starting later and/or ending earlier than the run |
| `split` | two or more |
| `not_reported` | none |

A Refined IBD segment that lies inside no run is listed as `refined_only`.
That would mean the two programs did not see the same markers, or that
Refined IBD does not work as described here.

The columns are the category, the run's samples, haplotypes, chromosome,
start, end, length, markers, carriers and frequency, then:
- the Refined IBD segments inside the run;
- every segment start or end inside the run, with the window edge it falls
  on, if any;
- the windows the run touches.

A summary of the counts goes to standard error.

## Tests

From the repository root:

```
python -m unittest discover -s tests
```

`tests/test_find.py` checks:
- the PBWT against a direct comparison of every pair of haplotypes on 300
  seeded random data sets;
- a hand-built case with every pair, including a sample's own two
  haplotypes, and its columns;
- the minimum length at its exact boundary and its default;
- runs on two chromosomes;
- the frequency format;
- the errors and the command line.

`tests/test_compare.py` gives each category once, with gzip and plain
inputs, and checks the window edges, the errors and the command line. All
data are synthetic. All 12 tests pass with Python 3.12 and 3.13.

## Licence

Every file in this repository except LICENSE is released under the GNU
General Public License, version 3 or (at your option) any later version
(SPDX: GPL-3.0-or-later). LICENSE is the text of that licence. The software
comes with no warranty (LICENSE, sections 15-16).

The modules `vcf`, `genetic_map`, `ids`, `input_it` and `java_compat` come
from [merge-ibs-segments](https://github.com/lakishadavid/merge-ibs-segments)
(commit `93402c9`). That is a modified Python translation of merge-ibd-segments
from Refined-IBD 17Jan20.102, Copyright (C) 2014-2016 Brian L. Browning,
GPL-3.0-or-later. Each module keeps his original notice and names the Java
files it was translated from. The BGZF reading in `input_it.py` is translated
from code Copyright (c) 2009 The Broad Institute under the MIT licence, whose
notice is kept in that file.

Copyright (C) 2026 LaKisha T. David.

## Not affiliated

This is an independent project. It is not written, maintained or endorsed by
Brian L. Browning, Richard Durbin or their institutions. Please report
problems in this repository's GitHub issues.

## Citation

If you use this software in published work, please cite this repository
with the version (tag or commit) you used, and the PBWT paper (Durbin 2014,
above). If you compare with Refined IBD, also cite Browning BL, Browning SR
(2013). Improving the accuracy and efficiency of identity-by-descent
detection in population data. *Genetics* 194(2):459-471.
doi:10.1534/genetics.113.150029. This is a request, not a condition of the
licence.
