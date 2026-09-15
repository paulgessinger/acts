# Sixteenth pass: reserve interval-bin storage once

Baseline: `d99b0fd900` (64 bins, 16 MiB scalar cache).

Counting bin memberships before filling the index avoids repeated vector growth.
The implementation computes each entry's first and last bins once, records those
ranges, reserves the exact bin capacities, and then appends original track
indices in their original order. Support bounds and all density arithmetic are
unchanged. Invalid intervals retain empty ranges.

## Performance

| Sample | Full-AMVF time reduction |
|---|---:|
| Initial 10-event A/B, five repetitions | 3.31% |
| Same 10 in the three-build run | 3.63% |
| Additional 90 events | 3.91% |
| All 100 events | **3.88%** |

In the 100-event comparison, the sum of per-event AMVF medians falls from
12797.624 to 12300.802 ms. Direct seeding changes by -0.71%.
The fresh profile's inclusive index-construction share falls from approximately
5.31% to 2.22% on event 3912208. Profile shares are overlapping samples, not
independently additive timing measurements.

The all-event run includes the original unoptimized implementation as a third
reference: its corresponding sum is 42327.263 ms, yielding **3.441× cumulative
full-AMVF speedup**. See the final report for the full configuration and memory
tradeoff. Process order rotates and reverses across events; every build has one
warmup and three measured repetitions, pinned to CPU 0. Compilation and Athena
reconstruction had finished before measurement.

## Validation and tooling

All six relevant suites pass: the four AMVF/density/Kalman suites plus
FullBilloirVertexFitter and IterativeVertexFinder, which also use TrackAtVertex.
All 100 events have exactly matching seed positions/widths and byte-identical
vertex/track dumps across all three builds. Run-3/Run-4 ODD, dense-z and
transformed-perigee stress comparisons also pass. The dense Run-3 stress timing
is 1.61% slower in this block; the performance improvement is workload-dependent.

The final naming of the bin-count constant leaves the compiled `.text` section
byte-identical to the measured candidate; its hashes are included in the data.
The replay driver now accepts an optional reference library/binary to reproduce
three-build comparisons. Both its existing two-build mode and new three-build
mode passed an exact-output smoke replay.

Raw results: [amvf_sixteenth_pass_results.json](amvf_sixteenth_pass_results.json).
