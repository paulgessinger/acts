# AMVF support-interval lookup optimization

Date: 2026-09-15

The Gaussian seeder now indexes the existing track support intervals into 64
z bins. On the ODD baseline workload this reduces complete finder/fitter runtime
by approximately half, with exactly matching reconstruction output.

## Controlled A/B results

| Athena profile | Baseline median (ms) | Optimized median (ms) | Speedup | Time saved |
|---|---:|---:|---:|---:|
| run4 | 185.52 | 96.18 | 1.93x | 48.2% |
| run3 | 155.22 | 84.66 | 1.83x | 45.5% |

The workload and configuration are those in [the baseline report](AMVF_ODD_BASELINE_REPORT.md):
1,258 ODD reconstructed tracks, 108 Run-4 vertices or 74 Run-3 vertices. No
selection, fitting, annealing, covariance, or seeding parameters were changed.

Measurements use the same benchmark executable with either the original or
optimized `libActsCore.so`, GCC 15.1.1, `RelWithDebInfo`, assertions disabled,
and affinity to logical CPU 0 on the Threadripper 3970X. There are three A/B
blocks per era, alternating which implementation runs first, each with two
warm-ups and eleven measured repetitions. The table reports the median of all
33 measurements per implementation. Serialization is outside the timed region.
Raw samples are in [amvf_optimization_timings.json](amvf_optimization_timings.json).

Run-4 baseline range: 184.75–186.52 ms; optimized: 94.51–97.00 ms.
Run-3 baseline range: 154.74–158.15 ms; optimized: 84.08–85.07 ms.

## Hotspot and implementation

The original Run-4 profile attributed 66.56% of cycles inclusively to seeding,
with 45.30% self time in `GaussianTrackDensity::addTrackToDensity`. Every trial
position scanned all remaining seeding tracks, including those outside their
existing support interval. Initially 538 tracks pass the seeder's selection;
only about 24 contribute at a typical track trial position in this fixture.

The new index stores track indices in every bin touched by their support
interval. It builds each bin's list in input order. A density query visits only
that bin's candidates, retaining the original strict interval check, Gaussian
expression, derivatives, summation order, trial positions, and maximum selection.
The bins select candidates; density is evaluated at the actual trial position.

Small inputs (fewer than 32 seeding entries), degenerate or non-finite trial
ranges, non-finite queries, and refinements outside the indexed range use the
original exhaustive scan. Infinite interval bounds are clamped into the index;
empty intervals and NaN bounds cannot contribute. Using the same monotonic
mapping during construction and lookup prevents gaps at rounded bin boundaries.
The fixed 64-bin limit bounds index storage and construction by O(nTracks).
Broad or densely overlapping intervals can limit the benefit.

The index is local to each `globalMaximumWithWidth` call and rebuilds from the
current entries, so state reuse and changing track lists cannot leave it stale.
The public state and configuration layouts are unchanged.

## Profile after optimization

A fresh Run-4 `perf record -F 499 --call-graph dwarf` comparison captured
4,038 baseline and 2,116 optimized samples, with no samples lost. Both used two
warm-ups and 41 measured finder calls. Inclusive shares overlap.

| Region | Baseline share | Optimized share |
|---|---:|---:|
| Seeding, inclusive | 67.58% | 37.35% |
| `addTrackToDensity`, self | 45.42% | 7.85% |
| Track-to-vertex compatibility association, inclusive | 11.88% | 23.01% |
| Fitter, inclusive | 10.93% | 22.62% |
| Fitter weight/update, inclusive | 7.74% | 16.53% |

Index construction accounts for about 2.1% inclusively in the optimized profile.
Scaling the seeding shares by the uninstrumented event medians suggests roughly
125 ms to 36 ms in seeding, or approximately 3.5x faster seeding. These subsystem
times are estimates from sampling; the event medians are directly measured.

Remaining self costs include trigonometry (`sincosf32x`, 9.72%), perigee position
conversion (`LineSurface::localToGlobal`, 6.54%), and parameter extraction/copies
(5.26%). Compatibility association and fitting each now occupy about 23% of the
shorter event. A next optimization pass can examine event-local caching of
track positions and avoiding repeated parameter copies. The fitter's remaining
cost is distributed across weight/update bookkeeping and small matrix work.

Local profiles and flat/inclusive reports are `baseline-perf.data`,
`final-perf.data`, `baseline-{flat,inclusive}.txt`, and
`final-{flat,inclusive}.txt` in the artifacts directory below.

## Correctness

- All three targeted suites pass: `TrackDensityVertexFinder`,
  `AdaptiveMultiVertexFinder`, and `AdaptiveMultiVertexFitter`.
- A new exhaustive-search reference checks exact maximum and width equality for
  Gaussian and parabolic refinement, randomized interval widths and track
  permutations, small and large lists, overlapping intervals, strict endpoints,
  bin boundaries, non-finite values, and refinement outside the lookup range.
- `--output-vertices` writes vertex positions, covariance, fit quality, track
  identities, weights, compatibility, fitted parameters, and fitted covariances
  at 17-digit precision. Baseline and optimized dumps match byte-for-byte in
  both eras for every A/B block.

## Additional workload checks

| Workload | Profile | Baseline (ms) | Optimized (ms) | Vertices |
|---|---|---:|---:|---:|
| shuffled | run4 | 189.372 | 98.504 | 108 |
| shuffled | run3 | 157.967 | 85.349 | 74 |
| small64 | run4 | 0.698 | 0.695 | 3 |
| small64 | run3 | 0.830 | 0.833 | 5 |
| dense | run4 | 4206.952 | 4215.146 | 8 |
| dense | run3 | 2763.643 | 2765.586 | 4 |

Each uses two warm-ups and seven measurements. All six baseline/optimized
output pairs match byte-for-byte. `shuffled` permutes the same event with Python
`random.Random(42)`; `small64` uses its first 64 tracks; `dense` scales all z0
values by 0.01, leaving the other parameters and covariance unchanged. The dense
case is an artificial algorithm stress test. It is not a physics event sample.
Its roughly 0.1–0.2% timing increase is below what this short check can resolve;
no speedup is claimed there. Likewise, the small input is essentially unchanged.
Raw samples are in [amvf_optimization_validation.json](amvf_optimization_validation.json).

## Reproduction

The isolated worktree is on `codex/amvf-optimization`, based on canonical ACTS
`d3f0323ac27ce7e5a97fccf9aa3be0c87787fb02` plus the baseline harness.
The original baseline worktree is preserved.

Local build paths:

- Original library: `/tmp/acts-amvf-build-main/lib64/libActsCore.so`
- Optimized build: `/tmp/acts-amvf-optimized-build`
- Input: `/tmp/odd-mu200-tracks.csv`
- Input SHA-256: `bc7274740f9afda8428cecae9c5f1ece0067996df0fd25826c4732d58ef81d5d`
- Detailed local artifacts: `/tmp/amvf-optimization-results`

For this pair of builds, whose public class layouts are unchanged:

```sh
LD_LIBRARY_PATH=/tmp/acts-amvf-build-main/lib64 taskset -c 0 \
  /tmp/acts-amvf-optimized-build/bin/ActsBenchmarkAdaptiveMultiVertexFinder \
  --input /tmp/odd-mu200-tracks.csv --athena-era run4 \
  --warmup 2 --repetitions 11 --output-vertices /tmp/baseline-vertices.txt

taskset -c 0 \
  /tmp/acts-amvf-optimized-build/bin/ActsBenchmarkAdaptiveMultiVertexFinder \
  --input /tmp/odd-mu200-tracks.csv --athena-era run4 \
  --warmup 2 --repetitions 11 --output-vertices /tmp/optimized-vertices.txt

cmp /tmp/baseline-vertices.txt /tmp/optimized-vertices.txt

ctest --test-dir /tmp/acts-amvf-optimized-build \
  -R '^(TrackDensityVertexFinder|AdaptiveMultiVertexFinder|AdaptiveMultiVertexFitter)$' \
  --output-on-failure
```

The optimized build uses the baseline's generated headers through
`ACTS_CODEGEN_PREBUILT_DIR=/tmp/amvf-prebuilt-codegen` to avoid network-dependent
code generation. Compiler and dependency locations are in its CMake cache.

## Scope

This is the existing controlled ODD particle-gun workload with 200 source
vertices and diagonal reconstructed covariance. It does not establish the
speedup for ATLAS production or ODD `ttbar + mu=200`. A representative multi-event
sample with full covariance remains the next physics/performance validation.

## Production representativeness audit

The measured speedup is conditional on the support-window overlap in the input.
The optimization adds no rejection cut: it indexes the bounds that the original
`GaussianTrackDensity` already checked. It optimizes the Gaussian seeder path;
other seeders do not use this index. The finder's `tracksMaxZinterval` (0.5 mm
in the audited Run-4 profile, 3 mm in Run 3) is a separate compatibility cut and
is unchanged.

For a symmetric d0/z0 covariance block with entries A = Cdd, B = Cdz,
D = Czz, completing the quadratic used by `addTracks` gives:

- Support centre: `z0 + d0 * B / A`.
- Half-width squared: `(D - B*B/A) * (4*S*S - d0*d0/A)`, where `S*S` is
  the effective `z0SignificanceCut`.

For diagonal covariance and small d0, the half-width is approximately
`2*S*sigma_z0`, hence about `24*sigma_z0` for the default S = 12. This factor
follows from the existing implementation's `constantTerm + 2*z0SignificanceCut`.
Correlations change both the centre and width; discarding them directly changes
which tracks overlap, and therefore changes the expected optimization benefit.

No beam-spot covariance enters that per-track support calculation.
`TrackDensityVertexFinder` passes only tracks to the density estimator, then
uses the constraint to construct the seed. Beam conditions still matter:
Athena selects tracks with a beam-spot reference, supplies the beam constraint
to the finder/fitter, and the luminous-region distribution controls the overlap
of track windows along z. Changing only the harness beam-spot sigma does not
reproduce a different collision distribution or detector resolution.

The local Athena checkout's `Projects/Athena/build_externals.sh` pins ACTS
47.5.0. Both `GaussianTrackDensity.hpp` and `GaussianTrackDensity.cpp` are
identical between that tag and our baseline ACTS main revision, so the original
support-window logic exists in the pinned version too. This does not establish
whole-chain equivalence or identify the release used by a particular campaign.

An additional configuration caveat: in ACTS 47.5.0 and this baseline, the
constructor computes `d0SignificanceCut` and `z0SignificanceCut` once. Athena's
wrapper subsequently assigns `d0MaxSignificance` and `z0MaxSignificance`, which
does not recompute those squared cuts. Defaults 3.5/12 agree; non-default
Gaussian property overrides require checking the effective cuts. A sensitivity
scan must initialize `Config(d0Sig, z0Sig)` or explicitly update both cuts.
No configuration behaviour was changed as part of this optimization.

Production validation remains open until the target campaign/release and input
are identified. Required validation is an A/B replay of selected production
tracks with their full covariance and perigee reference surfaces, matching
beam-spot position/covariance and geometry/field conditions, across multiple
Run-3/Run-4 events and occupancies. Record selected-track counts, support-window
width/overlap distributions, candidate fraction, fallback frequency, output
equality, and timing distributions. A realistic ODD ttbar + mu=200 sample is
useful for the ODD study, but cannot alone validate ATLAS detector resolutions
or production conditions.

## Suitable Athena ART inputs identified

Source inspected: local Athena checkout at `b4641de703a247b95e7d68ade2583a67a183c311`.
These are test definitions and accessible inputs; no reconstruction run or
production A/B validation has yet been performed.

Recommended Run-4 anchors in `InnerDetector/InDetValidation/InDetPhysValMonitoring/test/`:

- `test_run4_acts_vertex_PU200.sh`: dedicated vertex comparison, 100 events,
  `OnlyTrackingPreInclude` plus `actsLegacyWorkflowFlags` for its ACTS leg;
  produces `AOD.acts.root` and vertex IDPVM/DCube comparisons.
- `test_run4_acts_ttbar_PU200.sh`: use the ACTS leg with
  `OnlyTrackingRecoPreInclude,ActsConfig.ActsCIFlags.actsProductionFlags` for the
  production-oriented track population. This also enables ITk fast tracking.
  Its comparison plots focus on tracking, so retain the dedicated vertex test's
  vertex validation when using this configuration for the AMVF study.
- `test_run4_ttbar_PU200.sh`: standard Athena tracking with vertex validation,
  useful for an additional upstream-track population. Its old all-hadronic
  description disagrees with the currently resolved SingleLep input name.

All three resolve `defaultTestFiles.RDO_RUN4[0]` to:

```text
/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/PhaseIIUpgrade/RDO/ATLAS-P2-RUN4-04-00-00/mc21_14TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.recon.RDO.e8481_s4494_r16635/RDO.46493535.100evt.pool.root
```

The file resolves locally through CVMFS (578,528,447 bytes). The tests identify
it as PU200; the filename identifies a 14 TeV semileptonic ttbar sample and a
100-event subset. Event-level pile-up and beam conditions still need to be read
in the reconstruction job. The explicit conditions tag resolves to
`OFLCOND-MC21-SDR-RUN4-06`. The input directory identifies geometry
`ATLAS-P2-RUN4-04-00-00`; the tests let reconstruction read geometry from input
metadata.

For a sample with an explicit pile-up campaign configuration,
`Tools/TrfTestsART/test/test_trf_RUN4_r2a_ca_mt_mu200_overlay.sh` supplies 25
HITS-plus-presampled-background events and
`Campaigns.MC23PhaseIIPileUp200`. The campaign sets 200 collisions and a
190–210 pile-up profile. It is a heavier alternative because it includes overlay.

Run-3 companion: `test_ttbarPU40_reco.sh` in the same IDPVM test directory,
100 events from ART input dataset
`user.keli:user.keli.mc23a_13TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.recon.RDO.e8514_e8528_s4111_s4114_r14622_tid33359244_00`,
AMI tag `r14519`, run number `801271`, and current default conditions
`OFLCOND-MC23-SDR-RUN3-11-02`. Dataset availability has not been checked.
The generic `RDO_RUN3` default is low-mu minimum bias and should not replace
this ttbar input for the present performance study.

Next execution should derive its reconstruction command from the Run-4
production-flags ART leg, preserve the conditions setup, and capture selected
tracks with full covariance plus the perigee reference and actual beam-spot
conditions at the vertex-tool boundary. Compare baseline and optimized AMVF
on identical captured inputs before attributing a gain to production settings.

## Minimal track-parameter replay path

ACTS already has `CsvTrackParameterReader` in `Examples/Io/Csv`, which reads
per-event perigee parameters and the full 5x5 spatial covariance using
`TrackParameterData`. It initializes the unused time coordinate to zero and its
variance to one, matching the audited no-time Athena conversion. Its `beamspot`
configuration sets a common perigee reference point; it does not supply the
vertex constraint and is not event-dependent.

`RootAthenaNTupleReader` in `Examples/Io/Root` also reads the spatial covariance
and a beam-spot constraint, but expects a dedicated ntuple schema rather than
an arbitrary AOD. In its current implementation it assumes an origin perigee,
constructs only a diagonal beam-spot covariance, and leaves the time row/column
of its track covariance uninitialized. It should not be treated as a drop-in
production replay reader without validation/fixes. `RootTrackSummaryReader`
retains only diagonal track covariance and cannot close the current fidelity gap.
No exporter matching the Athena-ntuple reader's branch names was found in the
searched Athena tracking/ID/physics-analysis code or in the selected ART scripts.

The small export/replay interface should contain:

- Per event: event identity, mu, selected track order, actual beam-spot position
  and full 3x3 covariance, and the applicable field/configuration metadata.
- Per selected track: d0, z0, phi, theta, q/p; full 5x5 covariance; perigee
  reference point/transform and stable input identity.
- Explicit units, including the MeV-to-GeV q/p conversion and scaling of its
  covariance row and column.

Export after the same vertex-tool track selection and preferably immediately
after Athena's conversion to ACTS parameters. Exporting all AOD tracks would
otherwise change the input population. Retain a separate per-event metadata
record for the beam constraint and surface reference; the existing CSV reader
can supply the parameter representation but needs integration with that metadata.

The ART configurations produce AOD; their expensive upstream reconstruction is
needed only once if a matching AOD is not already available. Subsequent AMVF
measurements can run entirely from captured parameters. Hits, clusters and track
states are unnecessary for this replay. Truth and reference vertices are optional
validation data, not finder inputs.

The Gaussian seeder alone needs only track parameters/covariance and its
configuration, so its window-overlap speedup can be tested exactly without the
ATLAS field or tracking geometry. Complete AMVF also uses the field and a
propagator for compatibility and linearization. An isolated constant-field replay
can compare both implementations on realistic tracks, while exact Athena
whole-tool parity additionally requires matching that propagation setup.

## Executed ART export and replay (2026-09-15)

The proposed parameter-only replay is now implemented and exercised on 10
Run-4 PU200 ART events using the installed nightly's production configuration.
Full covariance export/import checks and exact A/B output comparisons pass for
every event. Aggregate Gaussian seeding is 2.70× faster; standalone full AMVF
is 1.73× faster. See [the ART results](AthenaExport/RESULTS.md) for exact sample
provenance, event timings, beam/window distributions, and the constant-field
limitation of the full-AMVF replay. This supersedes the earlier execution plan.
