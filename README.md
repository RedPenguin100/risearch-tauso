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
* Much faster (20x to 47x) by SIMD, caching and taking heavy calc out of loops
* Readability: functions that do one thing, named types, comments that explain why
* Tests: unit tests, end-to-end tests, and regression against a reference binary


**It is byte identical to the original, except the documented bugfixes**. You may read more in [HISTORY.md](HISTORY.md).

The main path that was optimized is what was previously called `RIs_linspace`. If the `force_start` path is of interest to you, please open an issue and similar changes can be performed.

## Performance

Against the original C, on an idle cluster node, one core, repeated on three nodes
which agree within 1%. Each timing is the best of three blocks, where a block runs
the binary N times and divides by N: one run alone is short enough that a single
`time` call would be measuring the 10 ms clock rather than the program.

| queries | targets | upstream C | 1.6.1 | speedup |
| --- | --- | --- | --- | --- |
| 16 x 20 nt | 4 000 records x 600 nt | 6.301 s | 0.263 s | **24.0x** |
| 16 x 20 nt | 960 records x 2 500 nt | 5.288 s | 0.144 s | **36.8x** |
| 16 x 20 nt | 120 records x 20 000 nt | 4.902 s | 0.109 s | **45.1x** |
| 16 x 20 nt | 1 record x 2 400 000 nt | 4.853 s | 0.131 s | **37.2x** |
| 64 x 20 nt | 4 000 records x 600 nt | 24.276 s | 1.125 s | **21.6x** |
| 64 x 20 nt | 960 records x 2 500 nt | 20.937 s | 0.578 s | **36.2x** |
| 64 x 20 nt | 120 records x 20 000 nt | 19.532 s | 0.415 s | **47.1x** |
| 64 x 20 nt | 1 record x 2 400 000 nt | 19.356 s | 0.437 s | **44.3x** |
| 256 x 20 nt | 4 000 records x 600 nt | 95.614 s | 4.616 s | **20.7x** |
| 256 x 20 nt | 960 records x 2 500 nt | 83.500 s | 2.311 s | **36.1x** |
| 256 x 20 nt | 120 records x 20 000 nt | 77.995 s | 1.634 s | **47.7x** |
| 256 x 20 nt | 1 record x 2 400 000 nt | 77.218 s | 1.660 s | **46.5x** |

Every target file holds the same 2 400 000 nt and differs only in how it is divided
into records, so the alignment work is the same in each. The 600 nt rows come out
lowest because a shorter record reports more hits, and every hit costs a traceback.

At 2 000 queries, where running the C would take about half an hour a row:

| queries | targets | 1.6.1 |
| --- | --- | --- |
| 2 000 x 20 nt | 4 000 records x 600 nt | 48.914 s |
| 2 000 x 20 nt | 960 records x 2 500 nt | 21.244 s |
| 2 000 x 20 nt | 120 records x 20 000 nt | 13.141 s |
| 2 000 x 20 nt | 1 record x 2 400 000 nt | 12.629 s |

Without AVX2 -- `RISEARCH_NO_AVX2=1`, which also turns off the batched sweep -- 16
queries against 960 x 2 500 nt takes 2.969 s, still 1.8x the C.

The inputs are generated from a fixed seed, so the table can be reproduced:

```
python3 RIsearch1/tests/performance/make_bench_data.py
RIsearch -q q_16.fa -t t_2500.fa -p3 > /dev/null
```

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
