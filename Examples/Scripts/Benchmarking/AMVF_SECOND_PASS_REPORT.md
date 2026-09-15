# AMVF second optimization pass

## Result

Relative to commit `b3b8082a23` (the support-interval index), this pass reduces
Gaussian seeding time by **11.94%** and standalone full-AMVF time by **5.54%**
on the same 10 Run-4 ART events. All 10 events improve in both measurements.

| Sum of per-event median times | Committed baseline | Second pass | Speedup |
|---|---:|---:|---:|
| One Gaussian seeding call | 49.218 ms | 43.341 ms | 1.136× |
| Full AMVF replay | 2158.226 ms | 2038.591 ms | 1.059× |

Each event uses one warmup and five measured repetitions, pinned to CPU 0.
Baseline/optimized order alternates across events. Both variants use the same
benchmark executable, compiler and configuration; only `libActsCore.so` changes.
The baseline library was preserved before rebuilding. Ratios use sums of
per-event medians. These are incremental gains over the first optimization.

## Profile and change

A fresh 199 Hz DWARF call-graph profile of ART event 3912208 (25 repetitions)
attributes about 70% inclusive time to Gaussian seeding, 12% to fitting and 10%
to compatible-track search. The per-candidate `addTrackToDensity` function
accounts for about 23% self time; exponential evaluation accounts for about
38% inclusive time. Percentages are approximate sampled costs, with nested
inclusive categories overlapping.

Move the private per-track density definition before the query loop and mark
it inline. This lets the compiler inline the candidate update and retain the
accumulators in the calling loop. Inspection of the rebuilt object confirms
that no `addTrackToDensity` call or standalone function remains. Gaussian
expressions, strict support bounds and track summation order are unchanged.
There are no new cutoffs, caches, approximations, configuration values, or
changes to public class layouts. The improvement is compiler/build dependent.

## Validation and scope

- Exact seed position/width equality on all 10 ART events.
- Byte-identical complete vertex/track dumps on all 10 ART events.
- Byte-identical Run-3 and Run-4 ODD fixture and dense-z stress outputs
  (three repetitions per variant/configuration).
- The TrackDensityVertexFinder, AdaptiveMultiVertexFinder and
  AdaptiveMultiVertexFitter test suites pass, including the exhaustive Gaussian
  oracle, small collections, interval edges and non-finite inputs.

The inputs, production selector, full covariance and actual beam constraints
are described in [the ART sample report](AthenaExport/RESULTS.md).
Full AMVF still uses the standalone constant 2 T field and simplified
propagation. This is not an in-Athena whole-tool timing measurement.

Raw timings, profile summaries and validation timings are retained in
[amvf_second_pass_results.json](amvf_second_pass_results.json).
Detailed perf captures and assembly are in `/tmp/amvf-pass2-results`.
