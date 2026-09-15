# Athena-default AMVF on an ODD Run-4-complexity workload

Date: 2026-09-15

Revision: `d3f0323ac27ce7e5a97fccf9aa3be0c87787fb02`, the fetched tip of
canonical `acts-project/acts` `main` at measurement time.

Athena reference: `154f86876b263a9e5b658535e37baa20bac6d643`, local
`atlas/athena` `origin/main` at configuration-audit time.

## Athena configuration audit

On current Athena `main`, ordinary offline Run-3 and Run-4 reconstruction
selects `VertexSetup.ActsGaussAMVF`. Specialized low-mu, minimum-bias,
beam-spot, vertex-lumi, high-pileup-pass, and heavy-ion modes select an
iterative finder and are not represented by this AMVF benchmark.

The harness now explicitly reproduces the following effective settings:

| Setting | Run 3 | Run 4 |
|---|---:|---:|
| Seeder | Gaussian track density | Gaussian track density |
| `tracksMaxZinterval` | 3 mm | 0.5 mm |
| Finder `minWeight` | 0.0001 | 0.02 |
| Finder `maxIterations` | 100 | 200 |
| Beam constraint | enabled | enabled |
| `doFullSplitting` | false | false |

Shared Athena component defaults are also set explicitly: annealing
temperatures `{1.0}` and cutoff 9; fitter iterations 30, relinearization
distance 0.5, relative-shift threshold 0.01, fitter minimum weight 0.001, and
smoothing enabled; Gaussian d0/z0 significance limits 3.5/12; finder track
significance 5, vertex chi2 18.42, merge significance 3, contamination 0.5,
real multivertex and fast compatibility enabled, no single-track vertices, no
seed constraint, no time fit, and initial variances `Vector4::Constant(1e8)`.
The impact-point estimator uses 20 iterations and precision `1e-10`.
In particular, the Run-4 finder weight of 0.02 is intentionally not copied to
the fitter: Athena leaves the fitter minimum weight at 0.001.

The audit trail in Athena is
`Tracking/TrkConfig/python/VertexFindingFlags.py` for setup and era-dependent
flags, `Tracking/Acts/ActsConfig/python/ActsPriVxFinderConfig.py` for effective
overrides, and
`Tracking/Acts/ActsVertexReconstruction/src/AdaptiveMultiPriVtxFinderTool.{h,cxx}`
for the component defaults and their mapping into ACTS.

Athena uses an `EigenStepper` with the ATLAS field and tracking geometry. The
harness now uses the same stepper implementation, but retains a constant 2 T
field and perigee-only geometry so it stays a small ODD benchmark. It uses the
Athena `BeamSpotCondAlg` fallback constraint at the origin with widths
`(0.15, 0.15, 53) mm`; production widths are conditions-dependent.

Track selection and Athena's xAOD conversion/final-output filtering are
outside the timed region. The CSV also retains only diagonal covariance terms,
whereas Athena passes the full covariance. These are known fidelity gaps; the
finder/fitter configuration itself is matched.

## Workload

The input is one ODD-reconstructed particle-gun event with 200 distinct
primary `z` positions and exactly four generated muons per position (800
generated particles). ODD digitisation and CKF produced 2,220 tracks; greedy
ambiguity resolution retained 1,251 tracks in the CSV summary and 1,258 fitted
parameter entries in `tracksummary_ambi.root`. The latter are replayed by the
benchmark; the Athena profiles reconstruct 108 vertices for Run 4 and 74
for Run 3.

This deliberately controlled event has the vertex multiplicity and track-count
scale needed for a Run-4 pile-up stress test. It is not yet the final physics
baseline: four uniform particle-gun muons per collision do not reproduce a
minimum-bias `ttbar + mu=200` track spectrum, vertex occupancy distribution, or
fake composition.

CSV parsing, parameter construction, finder construction, and the creation of
each fresh finder state are outside the timed region.

## Environment and timing

- AMD Ryzen Threadripper 3970X, pinned to logical CPU 0
- GCC 15.1.1, `RelWithDebInfo`, forced assertions disabled
- 2 warm-up runs followed by 11 measured runs
- Run 4 median: **184.02 ms/event**; mean 184.25 ms; range
  183.66--185.92 ms; 108 reconstructed vertices
- Run 3 median: **154.36 ms/event**; mean 154.30 ms; range
  153.86--154.70 ms; 74 reconstructed vertices
- Both runs replay 1,258 input tracks

An instrumented Run-4 repeat gave 187.85 ms median and the following
whole-process counters (13 finder calls, including warm-up): 10.18 billion
cycles, 14.35 billion instructions (1.41 IPC), 3.60 billion branches with
1.53% misses, and 279 million cache references with 8.37% misses. One-off
setup is negligible at this repetition count.

## CPU profile

`perf record -F 199 --call-graph dwarf` captured 496 cycle samples with none
lost. Inclusive percentages overlap and must not be added.

| Region | Inclusive | Self | Interpretation |
|---|---:|---:|---|
| `AdaptiveMultiVertexFinder::doSeeding` | 66.56% | <0.5% | Largest subsystem |
| `GaussianTrackDensity::addTrackToDensity` | 53.47% | 45.30% | Dominant leaf hotspot |
| `AdaptiveMultiVertexFitter::fit` | 11.04% | <0.5% | Fitter total |
| `AdaptiveMultiVertexFitter::setWeightsAndUpdate` | 7.69% | <0.5% | Main fitter subtree |
| `AdaptiveMultiVertexFinder::addCompatibleTracksToVertex` | 11.01% | 1.67% | Repeated compatibility scan |
| `AdaptiveMultiVertexFitter::doVertexSmoothing` | 2.72% | <0.5% | Smoothing contribution |
| `AdaptiveMultiVertexFitter::setAllVertexCompatibilities` | 0.63% | <0.5% | Compatibility bookkeeping |

Other visible self costs include vectorized `sincos` (4.80%),
`LineSurface::localToGlobal` (2.90%), copied track-parameter extraction through
the delegate (2.29%), and `exp` (2.06%).

## Initial conclusions

1. The first optimization target should be Gaussian seeding, not the fitter.
   `globalMaximumWithWidth` repeatedly evaluates every cached track entry at
   track-derived trial positions; its inner `addTrackToDensity` loop alone is
   45.3% of sampled cycles.
2. Under the Athena Run-4 profile, the fitter is only about 11% inclusive.
   Its time remains dispersed across weight/update bookkeeping, tree-map
   traversal, small dense linear algebra, and smoothing.
3. Compatibility association repeatedly computes track positions and
   trigonometry. Caching geometry-independent perigee quantities for the event
   is worth testing after seeding changes.
4. Before accepting an optimization, repeat this exact fixture for a clean
   A/B comparison and then validate it on an ODD `ttbar + mu=200` replay. The
   latter is required before generalising the absolute timing or hotspot mix to
   ATLAS Run 4.
