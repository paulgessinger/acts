# Run-4 ART parameter replay: measured results

The interval index improves Gaussian seeding on tracks reconstructed with the
Run-4 PU200 ART production configuration. Across 10 events, the sum of per-event
median times falls from **131.11 to 48.50 ms (2.70×)** for one Gaussian seeding
call, and **3712.93 to 2140.74 ms (1.73×)** for standalone full AMVF.
All seed positions/widths compare exactly; all fitted vertex/track dumps are
byte-identical between baseline and optimized libraries.

## Sample and configuration

Executed 2026-09-15 with Athena 25.0.73 nightly `2026-09-14T2100`, GCC 15.2,
and the production-flags leg of
`InDetPhysValMonitoring/test/test_run4_acts_ttbar_PU200.sh`:
`OnlyTrackingRecoPreInclude,ActsConfig.ActsCIFlags.actsProductionFlags`.
This installed nightly uses ACTS 47.7.0; the A/B replay uses the two ACTS builds
described in the parent optimization report. Their public layouts are unchanged.

The installed release's TestDefaults resolves geometry `ATLAS-P2-RUN4-05-00-00`
and conditions `COND-MC21-SDR-RUN4-06`. This differs from the older local Athena
checkout defaults audited earlier. The full resolved input RDO path, release,
and conditions are retained in [art_run4_results.json](art_run4_results.json).
Reco_tf completed with exit code 0 and validated its 10-event AOD.

All 10 events have actual and average mu=192. The production vertex selector
retains 907–1471 tracks per event (12,332 total), with no missing covariances.
The beam position is (0,0,0) and its standard deviations are (0.012,0.012,50) mm.
All selected perigee reference transforms are the identity. The export keeps the
full correlated 5×5 covariance, original track order, and q/p unit conversion.
All six reconstructed parameters and all 36 covariance elements compare exactly
against the exported representation, including the initialized time coordinate.

The Gaussian's existing d0 significance cut leaves 863–1382 seed tracks per
event. With its default z0 significance configuration (12), per-event median
support half-widths are 5.27–6.54 mm, and 95th percentiles are 55.54–68.29 mm;
initial z positions span 235–289 mm. Thus realistic covariance tails do widen
some windows, but enough windows remain local for the index to reduce scanning.
No Gaussian cutoff or production selection cut was tightened.

## Event results

Milliseconds, median of five measured repetitions after one warmup, pinned to
CPU 0. Baseline/optimized process order alternates across events; input parsing
and output serialization are outside timing. Each pair uses the same executable
with only the linked ACTS library changed. Aggregate speedups above are ratios
of sums of per-event medians, not averages of event speedup ratios.

| Event | Selected tracks | Seeder baseline → optimized (ms) | AMVF baseline → optimized (ms) | Vertices |
|---|---:|---:|---:|---:|
| 3912195 | 1073 | 9.10 → 3.15 | 250.82 → 150.19 | 86 |
| 3912203 | 907 | 6.63 → 2.81 | 183.08 → 119.15 | 82 |
| 3912208 | 1447 | 17.58 → 6.20 | 544.29 → 292.98 | 117 |
| 3912201 | 1431 | 16.85 → 6.06 | 533.98 → 279.77 | 116 |
| 3912211 | 1118 | 9.79 → 3.62 | 292.54 → 177.30 | 80 |
| 3912207 | 1177 | 11.91 → 4.39 | 295.86 → 177.88 | 91 |
| 3912204 | 1331 | 15.56 → 5.56 | 436.77 → 255.40 | 98 |
| 3912206 | 1471 | 19.59 → 7.94 | 520.40 → 308.70 | 105 |
| 3912205 | 1044 | 9.24 → 3.37 | 225.53 → 141.38 | 86 |
| 3912209 | 1333 | 14.86 → 5.40 | 429.65 → 237.99 | 109 |

## Scope and reproducibility

The seeder needs no field or detector geometry, so these measurements directly
test the window optimization with the ART-selected parameters and covariances.
Full AMVF uses the actual exported beam constraint and Run-4 finder settings,
but retains the benchmark's constant 2 T field and simplified propagation.
Its timing and exact A/B agreement are **standalone replay results**, not an
in-Athena whole-tool speedup or a comparison to Athena's reconstructed vertices.
Ten events from one Run-4 MC setup do not establish performance for every
production campaign; Run-3 and other beam/track configurations remain unmeasured.

The original AOD, event CSVs/metadata, effective flags, configuration pickle and
transform logs are retained at
`/home/pagessin/dev/acts/projects/amvf-art-data/run4-production-10`.
Reproduction commands and the export contract are in [README.md](README.md).
The three targeted ACTS vertexing test suites also pass after adding the replay
importer. Synthetic correlated-covariance and translated-reference input passes
roundtrip verification and both A/B output checks.
