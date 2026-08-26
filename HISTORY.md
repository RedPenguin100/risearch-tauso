# Change history from Original C

Walkthrough of all changes in this repository from the original C

## Focus

We focused on the performance RIs_linSpace variant of the code. The most common use case tested was
for short DNA strands attaching to RNA strands (to mimick an ASO).

The force_start path was slightly improved but not routinely tested against for performance.
In the future, we might improve that as well. (Please open issues if that is of your concern)

## Benchmarking

The code now is about 30x faster than the original C in various settings: 

| queries | target file | upstream C | 1.6.0 | speedup |
| --- | --- | --- | --- | --- |
| 16 | 4 000 records x 600 nt | 5.100 s | 0.175 s | **29.0x** |
| 16 | 1 000 records x 2 500 nt | 4.745 s | 0.150 s | **31.6x** |
| 16 | 120 records x 20 000 nt | 4.321 s | 0.134 s | **32.2x** |
| 16 | 1 record x 2 400 000 nt | 4.241 s | 0.165 s | **25.6x** |

the main speedups were in the traceback and the linSpace algorithms, 
and for multi-target setting some caching is often needed to retrieve this speed.

The optimal speed is gained per query when sending queries as multiples of 16.
Most of the improvements rely on AVX2, but not all. The scalar code is about 1.4x faster
even on machines without AVX2.

## Backwards compatibility

For the 1.2.0 to 1.6.0 versions the output of this program is mostly byte-identical 
to the original. Where it deviates, is documented and those are specific bug-fixes that 
we addressed [BUGFIXES.md](BUGFIXES.md).

---
## Before the fork

RIsearch1 as written by Anne Wenzel and extended by Giulia Corsi: one C file,
three DP matrices, a linear-space sweep to find where alignments end, and a
full-matrix traceback to recover the pairing.

## 1.2.0 — 2026-05-28 - python package

Created a python package `risearch-tauso` and nothing else

## 1.3.0 — 2026-08-10

* Instead of accessing DSM matrix directly, QueryProfile was created for better cache usage in hotpath and performance boost
* Replaced slow `isalpha` with our implementation for a small gain
* Began moving the code to C++, fixed memory leaks, const correctness etc
* Compile options - changed visibility to hidden, added LTO and gc-sections
* Removed the possibility of 0 length record to prevent bugs


## 1.4.0 — 2026-08-11

* Began integrating AVX2 into the project, for a significant performance boost, cumulatively 3x faster than original C. Integration was for faster single-query iteration
* Fixed a bug where the header would re-read the last record in an infinite loop

## 1.4.1 – 1.4.3 — 2026-08-11 / 2026-08-12

* Swapped traceback loop order for quicker iteration for a minor performance boost, consistent looping with main loop

## 1.5.0 — 2026-08-13

* Added AVX2 to traceback
* Moved many AVX2 primitives to `avx2/primitives.h` for making it easier to read the code
* Switched printf to the faster and safer fmt library (same output)
* Added `has_positive_gap` which reverts runs to scalar for backwards compatibility
* Sequence encoding was moved from switch to lookup table also for a small performance gain

## 1.6.0b1 — 2026-08-21

* Added possibility of using int16 instead of int32 in the main matrices (M, Ix, Iy) for better SIMD
* "BatchedQuery" introduced - running AVX2 on multiple queries at once instead of a single query. Big performance boost

## 1.6.0b2 — 2026-08-22

Older compilers (gcc10) would default to swapping 256-bit access into two 128-bit halves, which significantly slowed 
the code. Later compilers (gcc11) changed that, so the flag only matters for older toolchains

* Added a compile options to remove deprecated default from GCC10, for an about 45% performance boost(!!!)

## 1.6.0 - 2026-08-26

Readability and correctness:
* `getopt` result was held in char, we switched to int
* Fixed a `realloc` leak in ByteBuffer
* Removed the now dead code regarding matrices
* Began building CI with GCC10 as well in manylinux2014 image (wheel publishing target)

Performance improvements:
* Cache the QueryProfile, the BatchedQueryProfile
