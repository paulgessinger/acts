# Eighth pass: reuse track weights across Kalman updates

Baseline: `689dee41cd`.

## Results

On the 10 Run-4 ART events, the sum of per-event full-AMVF medians falls from
1370.315 to 1289.490 ms: **5.90% less time (1.063× speedup)**. Seeding time is
unchanged within noise (43.903 versus 43.791 ms). All seed results and full
vertex/track dumps match exactly.

An initial inline-cache prototype improved time by 5.16%, but raised peak RSS
on event 3912208 from 22,180 to 26,404 KiB. The final version allocates a cache
only on first fitting use and lets track copies share unchanged cache storage.
Its measured peak RSS is 23,248 KiB, about 1 MiB above baseline. These are
single process-level measurements, not universal memory bounds.

## Change and correctness

TrackAtVertex retains the covariance block inverse used by the Kalman updater.
The cache compares all covariance coefficients by their bit representations
and records whether five or six parameters were inverted. Any covariance or
dimension change forces recomputation using the original Eigen expression.
The updater reuses the result across annealing iterations and during smoothing.
This avoids relying solely on a relinearization flag or on caller notifications.

TrackAtVertex copies share a read-only cache; an update detaches before changing
a shared cache. Cache access is non-const, in the same exclusive track-update
context that already mutates TrackAtVertex. There is no global shared cache.

New tests exercise repeated requests, correlated time covariance, switching
between 3D and 4D fits, covariance mutation, copied-track independence and
replacement of the complete linearized state. All four targeted vertexing
suites pass, including the 4D fitter and smoothing test.

## Measurement

One warmup and five samples per ART event, CPU 0, alternating A/B order. Since
TrackAtVertex's layout changes, each variant uses its own executable built
against its own headers and library. The replay driver's baseline-binary option
keeps this comparison ABI-safe. This remains standalone constant-2-T replay,
with the same ART input and beam conditions described in the earlier reports.

Raw data: [amvf_eighth_pass_results.json](amvf_eighth_pass_results.json).
Logs, memory observations and profile: `/tmp/amvf-pass8-results`.

Final shared-cache Run-3/Run-4 ODD and dense-z stress replays also have
byte-identical outputs (three repetitions per variant/configuration).
