# AMVF third pass: reuse unchanged Gaussian sums

Baseline: `85a18896ec` (the committed inlining optimization).

## Results

Two complete timing blocks compare the cache with the committed second pass.
Numbers are sums of per-event medians across the same 10 ART events.

| Block | Measurement | Baseline (ms) | Cached implementation (ms) | Time reduction |
|---|---|---:|---:|---:|
| 1 | seeder | 43.476 | 45.760 | -5.25% |
| 1 | amvf | 2036.630 | 1637.787 | 19.58% |
| 2 | seeder | 43.403 | 44.100 | -1.61% |
| 2 | amvf | 2046.551 | 1645.149 | 19.61% |

Full-AMVF improves on every event in both blocks. The direct single-call
seeder has a small overhead in the repeat block; it does not enable caching.
Block 1 includes an anomalously slow optimized seeder measurement on event
3912207 (6.10 versus 4.06 ms), which motivated the complete repeat block.
Both blocks are retained rather than discarding the outlier.

Raw timings and stress results are in
[amvf_third_pass_results.json](amvf_third_pass_results.json).

## Mechanism

The Gaussian density is repeatedly evaluated while AMVF removes seed tracks.
A temporary diagnostic on ART event 3912208 observed 83,089 repeated query / 
contributor-signature pairs among the first 150,000 queries (55.4%), representing
47.5% of evaluated contributions. This diagnostic used a hash to estimate reuse;
the implementation uses exact comparisons, not hashes. Instrumented runtimes
were excluded from performance results.

An opt-in cache stores the full density and its two derivatives at finite,
nonzero z positions. Before the next search it compares all six double fields
of each density entry by their bit representations. If the new entries are an
ordered subsequence of the previous entries, it erases cached queries strictly
inside every removed entry's support interval. All remaining query sums have
exactly the same contributing sequence. Additions, reordering, or coefficient
changes reset the cache. Affected sums are evaluated from scratch in input
order; no contributions are subtracted from previous sums.

The cache holds at most eight queries per current density entry, clearing when
shrinking input exceeds that bound. Copied states detach before mutation.
Non-finite queries and signed zero bypass the cache. Cache storage is local to
the caller-owned state; it is not global or shared between concurrent events.
TrackDensityVertexFinder keeps this state across searches and rebuilds the
coefficients from the actual supplied tracks each call, avoiding assumptions
about track-pointer lifetime or removal notifications.

Direct GaussianTrackDensity users retain the uncached default. Persistent
vertex-finder state enables caching explicitly; a cold search pays cache
construction costs to benefit subsequent searches.

## Validation and measurement

The public State layouts change, so each A/B variant uses its own executable
built against the matching headers and library. The replay driver now accepts
`--baseline-binary` for this purpose. Both builds use the same compiler and
configuration. ART timings use CPU 0, one warmup and five measured repetitions
per event, with alternating baseline/optimized order.

The ART input sample and the constant-field limitation of full-AMVF replay
remain as described in [the sample report](AthenaExport/RESULTS.md). These are
standalone replay measurements, not an in-Athena whole-tool benchmark.

## Correctness and stress coverage

All 10 ART events have exact seed position/width equality and byte-identical
full vertex/track dumps. The three targeted vertexing suites pass. New oracle
coverage exercises repeated calls, removals, copied-state isolation, reordering,
changed coefficients, additions, empty collections and state reuse, for both
Gaussian and parabolic refinement.

Run-3 and Run-4 ODD and dense-z stress replays also have byte-identical outputs
(three repetitions per variant/configuration). ODD full-AMVF times improve
6.2% and 13.7%, respectively. Dense-z times change by less than 1%; removing
tracks with broadly overlapping windows leaves little reusable work.
