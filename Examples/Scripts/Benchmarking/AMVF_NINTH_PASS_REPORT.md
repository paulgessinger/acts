# Ninth pass: bounded reuse of individual Gaussian values

Baseline: `57fb2c72c2`.

## Results and memory tradeoff

The 10-event ART comparison reduces aggregate full-AMVF time from 1291.969 to
1218.786 ms: **5.66% less time (1.060×)**. Direct seeding improves 1.10%
(43.793 to 43.312 ms). All seed results and full vertex/track dumps match exactly.

Peak RSS on ART event 3912208 is 23,040 KiB baseline and 31,928 KiB optimized:
**about 8.7 MiB additional process memory**. Scalar table storage is capped at
8 MiB per density state, plus row/index/container overhead. Small events use
smaller tables. This round trades bounded memory for CPU time; the bound is not
a claim that total process memory increases by at most 8 MiB.

## Mechanism

Full-query sums must be invalidated when a contributing track is removed, but
the surviving tracks' individual exponentials at the same z are unchanged.
The cache retains scalar exponential values at selected initial trial z values,
indexed by their original density-entry slots. It uses at most 1024 rows and
at most 1,048,576 doubles in total. Rows are allocated and populated on demand.

Ordered-subsequence matching maps surviving entries to original slots. It uses
the existing bitwise comparison of all six entry fields. Additions, coefficient
changes or reordering discard the table; removals preserve surviving values.
Copied density states copy the table when detaching. Non-finite and signed-zero
queries bypass it. Cached exponentials still pass the original strict support
check and are accumulated in original track order, with unchanged derivative
arithmetic. There is no subtraction, approximate exponential or reassociation.

An initial shared-loop prototype slowed uncached seeding. The final version
keeps a separate compact uncached evaluation loop, eliminating that regression.

## Validation and measurement

Four targeted vertexing suites pass. The exhaustive density oracle now runs
with caching both disabled and enabled, including randomized widths/orders,
small/large collections, interval boundaries, non-finite values and refinements
outside the indexed range. Cache mutation and copied-state tests also pass.
Run-3/Run-4 ODD and dense-z replays have byte-identical outputs.

ART timings use CPU 0, one warmup and five samples per event, alternating A/B
order. Public object layouts are unchanged, so the same executable loads each
matching library. The same ART inputs and constant-2-T standalone propagation
limitation apply. Raw results:
[amvf_ninth_pass_results.json](amvf_ninth_pass_results.json).
