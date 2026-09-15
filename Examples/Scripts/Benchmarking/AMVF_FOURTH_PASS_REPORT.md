# AMVF fourth pass: reuse the Kalman track weight matrix

Baseline: `b0339da0c6` (the committed Gaussian query cache).

## Result

Across the same 10 Run-4 ART events, the sum of per-event full-AMVF median times
falls from **1638.730 to 1576.773 ms**, a **3.78% time reduction (1.039×)**.
Every event improves. Gaussian seeding is unchanged within timing noise
(43.871 versus 43.891 ms aggregated). All seed results compare exactly and all
full vertex/track dumps are byte-identical.

Each event uses one warmup and five measured repetitions, pinned to CPU 0.
The replay driver alternates baseline/optimized order across events. Both
variants use the same executable and corresponding ACTS shared libraries.
The change only enlarges a function-local detail cache, leaving the benchmark's
public state layouts unchanged. Baseline libraries were saved before rebuilding.

## Change

`KalmanVertexUpdater::detail::calculateUpdate` and `trackParametersChi2` run
successively on the same linearized track. Both previously inverted the same
covariance block: 5×5 for a spatial vertex and 6×6 for a spacetime vertex.
The updater now stores that inverse in its existing per-update cache and reuses
it for the track chi2 calculation. The inversion expression and downstream
matrix expressions remain unchanged.

The inverse is recomputed for every update, so reuse needs no invalidation
logic across tracks, annealing iterations or relinearization. It uses the
appropriate covariance block directly, retaining correct treatment of time
correlations instead of taking a block of the existing full inverse.

## Validation and scope

The three targeted vertexing suites pass, including AdaptiveMultiVertexFitter's
4D time-fitting and smoothing test. ART exact-output comparisons cover all
10 events. The previous profile attributed about 8% of total sampled time in
ART event 3912208 to the 5×5 inversion routine; this is an inclusive sampled
cost, not a prediction that all of it can be removed by this change.

The ART provenance and standalone constant-2-T propagation limitation remain
as described in [the sample report](AthenaExport/RESULTS.md). These results
measure isolated replay, not an in-Athena whole-tool speedup.

Raw timings and validation measurements are retained in
[amvf_fourth_pass_results.json](amvf_fourth_pass_results.json).
Build/test logs, output dumps and perf capture: `/tmp/amvf-pass4-results`.

The dedicated KalmanVertexUpdater suite also passes. Run-3 and Run-4 ODD and
dense-z stress replays (three repetitions per variant/configuration) have
byte-identical vertex/track outputs.

The follow-up profile on ART event 3912208 reduces the 5×5 inversion routine's
inclusive share from 7.97% to 5.51%, consistent with removing the duplicated
inversion. Dense-z stress timings improve by roughly 18–19%; those artificial
inputs spend a larger fraction of time in fitting than the ART sample.
