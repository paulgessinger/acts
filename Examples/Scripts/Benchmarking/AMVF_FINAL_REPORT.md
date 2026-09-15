# AMVF optimization: final 100-event ART validation

The final implementation (`e38c13787c`) is **3.44× faster for full
AMVF** and **3.09× faster for one-off Gaussian seeding** than the
original core. Full-AMVF runtime falls **70.9%**. All
100 events have exactly matching seed positions/widths and byte-identical fitted
vertex/track output across the original, previous optimized and final builds.

## Direct cumulative measurement

Times below are sums of 100 per-event medians in milliseconds. Speedups are
ratios of those sums, not averages of event-by-event ratios or products of the
individual optimization-round gains.

| Workload | Original core | Previous optimized | Final | Original / final |
|---|---:|---:|---:|---:|
| Full AMVF | 42327.263 | 12797.624 | 12300.802 | 3.441× |
| One-off seeding | 1420.703 | 463.807 | 460.509 | 3.085× |

The additional 90 events show **3.45× full-AMVF speedup**. All 100
also participated in the final cache-budget validation; these cover one ART
sample and are not an independent production campaign.

The original library is from `d3f0323ac27ce7e5a97fccf9aa3be0c87787fb02`.
Its compatible benchmark frontend is from `85a18896ec`, which supplies the ART
CSV/metadata reader while retaining the original public layouts. The previous
optimized library is `d99b0fd900`. Each executable uses its matching library;
public-layout changes are never crossed by substituting an incompatible library.

GCC 15.1.1, RelWithDebInfo, assertions disabled, Threadripper 3970X, CPU 0.
Each build has one warmup and three measured repetitions per event. Three-build
execution order rotates and reverses across events. Reconstruction and compilation
finished before timing; parsing and output serialization are outside the timed
region. Raw samples are in the `matrix` field of
[amvf_sixteenth_pass_results.json](amvf_sixteenth_pass_results.json).

## Production-derived inputs

The complete 100-event Run-4 ttbar PU200 ART RDO was reconstructed with Athena
25.0.73 nightly `2026-09-14T2100`, using the production-flags leg of
`test_run4_acts_ttbar_PU200.sh`: `OnlyTrackingRecoPreInclude` and
`actsProductionFlags`. Reco_tf exited successfully and validated its 100-event
AOD. Geometry is `ATLAS-P2-RUN4-05-00-00`; conditions are
`COND-MC21-SDR-RUN4-06`.

The production vertex selector retained **128,984 tracks**, 891–1,738 per event,
with actual pileup 191–209 and no missing covariances. Beam standard deviations
are approximately (0.012, 0.012, 50) mm. The Run-4 replay uses the exported beam
constraint, 0.5 mm track z interval, 0.02 finder minimum weight and 200 maximum
finder iterations. No Gaussian significance cut, selection, annealing or fitting
parameter was tightened to obtain these gains.

Every exported parameter and all 36 entries of each reconstructed covariance
matrix passed the importer check, including the initialized time entries. The
first ten track/metadata pairs are byte-identical to the earlier export. Input
hashes, event metadata, effective flag records and configuration-file hashes are
in [amvf_final_validation.json](amvf_final_validation.json). The AOD, CSVs,
metadata, flags, configuration pickle and logs remain under
`/home/pagessin/dev/acts/projects/amvf-art-data/run4-production-100`.

**Scope:** full AMVF still uses the standalone constant 2 T field and propagation
setup. These are exact A/B replay comparisons on ART-derived inputs. In-Athena
whole-tool timing and equality to Athena's reconstructed vertices are unmeasured.

## Retained changes and memory

- Index existing Gaussian support intervals into 64 bins, with inline accumulation.
- Reuse unchanged complete density queries and bounded individual exponentials;
  preserve strict bounds, original summation order and exact invalidation.
- Pool query-cache allocations and reserve interval-bin storage before filling.
- Compute global input z positions once per find, after the first usable seed.
- Reuse covariance inverses and track-only Kalman projection matrices, with
  covariance/Jacobian/dimension checks and independent mutation of copied caches.

The scalar exponential table is capped at **16 MiB per density state**, plus
container overhead. The 32 MiB experiment was only 0.22% faster across 100 events,
so the smaller limit was retained. Peak process RSS, measured separately with
one warmup and three repetitions:

| Event | Selected tracks | Original RSS | Final RSS |
|---|---:|---:|---:|
| 3912208 | 1,447 | 21.18 MiB | 38.00 MiB |
| 3912243 (largest track count) | 1,738 | 23.04 MiB | 42.50 MiB |

These are whole-process peaks for those two events, not a bound on process or
Athena memory. The 16 MiB limit applies only to retained scalar table entries.

## Verification and remaining profile

All six relevant suites pass: AdaptiveMultiVertexFinder,
AdaptiveMultiVertexFitter, TrackDensityVertexFinder, KalmanVertexUpdater,
FullBilloirVertexFitter and IterativeVertexFinder. Coverage includes exhaustive
Gaussian comparisons, cache invalidation/copies, 3D/4D fitting and time
correlations. Run-3/Run-4 ODD, dense-z and transformed-perigee stress comparisons
also retain exact output. Runtime gains are workload-dependent; the final round's
dense Run-3 stress timing was 1.61% slower in its three-sample block.

The final 499 Hz DWARF profile on event 3912208 attributes about
68.5% inclusively to Gaussian seeding,
21.9% to scalar exp and
2.2% to index construction
(down from 5.3% before reserving bins). Inclusive shares overlap. The full report
is [amvf_final_profile.txt](amvf_final_profile.txt).

Rejected experiments are recorded in rounds 6, 10, 12, 13 and 15: insertion hints,
refinement helper inlining, exact-order derivative batching, different bin counts
and the full-bin contiguous-scan shortcut. Remaining time is dominated by scalar
Gaussian/refinement arithmetic and fitting. No further clearly beneficial,
low-risk exact-output change was identified in this review. Approximate math or
changed seeding/selection semantics were outside this optimization's scope.

## Reproduce

Use [the export/replay instructions](AthenaExport/README.md). For the final
three-build comparison, set the baseline to the `d99b0fd900` library and add the
original reference to the final executable's replay:

```sh
python3 Examples/Scripts/Benchmarking/AthenaExport/replay_exports.py \
  /tmp/acts-amvf-optimized-build/bin/ActsBenchmarkAdaptiveMultiVertexFinder \
  /home/pagessin/dev/acts/projects/amvf-art-data/run4-production-100/tracks \
  --baseline-library /tmp/amvf-pass16-baseline/lib64 \
  --optimized-library /tmp/acts-amvf-optimized-build/lib64 \
  --reference-binary /tmp/amvf-pass3-baseline/ActsBenchmarkAdaptiveMultiVertexFinder \
  --reference-library /tmp/acts-amvf-build-main/lib64 \
  --repetitions 3 --output /tmp/amvf-final-replay.json
```
