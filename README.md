# RIsearch1 (tauso fork)

RIsearch1: RNA–RNA, RNA-DNA and DNA-DNA interaction prediction using a simplified
nearest-neighbor energy model.

This fork carries **RIsearch1 only**. RIsearch2 and the siRNA off-target pipeline
live upstream at [RTH-tools/risearch](https://github.com/RTH-tools/risearch); they
were removed here because nothing downstream of this fork uses them. Their history
is still in this repository if you need it.

## What this fork changes

The fork is actively maintained and covers the RIsearch1 executable, above all
the standard linear-space search path. Its goals are:

* Correctness: fix the crashes, leaks and hangs listed in [BUGFIXES.md](BUGFIXES.md)
* Much faster (30x) by SIMD, caching and taking heavy calc out of loops
* Readability: functions that do one thing, named types, comments that explain why
* Tests: unit tests, end-to-end tests, and regression against a reference binary


**It is byte identical to the original, except the documented bugfixes**. You may read more in [HISTORY.md](HISTORY.md).

The main path that was optimized is what was previously called `RIs_linspace`. If the `force_start` path is of interest to you, please open an issue and similar changes can be performed.

## Performance

Against the original C - on an idle cluster node, one core, best of five, repeated on three nodes:

| queries | targets | upstream C | 1.6.0 | speedup |
| --- | --- | --- | --- | --- |
| 16 x 20 nt | 4 000 records x 600 nt | 5.100 s | 0.175 s | **29.0x** |
| 16 x 20 nt | 1 000 records x 2 500 nt | 4.745 s | 0.150 s | **31.6x** |
| 16 x 20 nt | 120 records x 20 000 nt | 4.321 s | 0.134 s | **32.2x** |
| 64 x 20 nt | 1 000 records x 2 500 nt | 19.503 s | 0.624 s | **31.2x** |

Every target file holds the same ~2.4 Mb of sequence and differs only in how it
is divided into records, so the alignment work is the same in each.


## Python package

A precompiled binary from this fork is published to PyPI as `risearch-tauso`, so
RIsearch is a normal Python dependency for tauso and any other downstream tool:

```bash
pip install risearch-tauso
```

```python
import risearch_tauso, subprocess
subprocess.run([risearch_tauso.executable_path(), "-q", "query.fa", "-t", "target.fa"])
```

Or as a CLI shim, which forwards all arguments straight to the bundled binary:

```bash
risearch-tauso -q query.fa -t target.fa
python -m risearch_tauso -q query.fa -t target.fa
```

The PyPI package is **not** canonical upstream RIsearch — it is the tauso-team
fork. Use upstream if you want the unmodified tool.

## Differences from upstream

Bugs fixed relative to upstream are recorded in [BUGFIXES.md](BUGFIXES.md), and
the release-by-release story is in [HISTORY.md](HISTORY.md).

## Building from source

Requires a C++17 compiler and CMake 3.20 or newer.

```bash
cmake -S RIsearch1 -B RIsearch1/build -DCMAKE_BUILD_TYPE=Release
cmake --build RIsearch1/build -j
```

That produces `RIsearch1/bin/RIsearch` and `RIsearch1/bin/RIsearch.dbg`, the
second with the debug tracing compiled in. The wheel build runs the same CMake via
`setup.py`.

## Tests

```bash
./RIsearch1/build/tests/risearch1_tests
```

Unit tests for the nucleotide coding, the max helpers, the alignment symbols, the
energy matrix, the int16 bound, the query profiles, the batched sweep and the
FASTA reader, including death tests for the inputs that are refused; end-to-end
tests that run `main()` in process with a constructed argv and assert on its
output; and throughput benchmarks under `Performance.*`. Pass
`-DRISEARCH_BUILD_TESTS=OFF` to skip googletest entirely, as the wheel build does.

`RISEARCH_NO_AVX2=1` in the environment forces the scalar kernels, which is how
the vector and scalar paths are compared from a single build.

## Running

```bash
RIsearch -q query.fa -t target.fa
```

Both files may hold several sequences; RIsearch scans all against all. Single
sequences can be given directly with `-Q acgu -T acgu`. See `RIsearch1/Manual.pdf`
for the full option list.

## Copyright

Copyright 2021 by the contributors; see `RIsearch1/README`.

RIsearch1 is released under the GNU General Public License version 3. This is free
software: you can redistribute it and/or modify it under the terms of that licence,
either version 3 or (at your option) any later version. You should have received a
copy of the GNU General Public License along with RIsearch — see the file COPYING.
If not, see <http://www.gnu.org/licenses/>.

This software is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citation

If you use RIsearch in a publication, please cite:

**RIsearch: fast RNA-RNA interaction search using a simplified nearest-neighbor
energy model.** Wenzel A, Akbasli E, Gorodkin J. *Bioinformatics*. 2012 Nov
1;28(21):2738-46.

## Contact

For problems with this fork, open an issue here. For upstream RIsearch:
<software+crispron@rth.dk>
