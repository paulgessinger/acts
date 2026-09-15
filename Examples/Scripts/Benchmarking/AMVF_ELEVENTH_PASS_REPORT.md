# Eleventh pass: cache track-only Kalman projection matrices

Baseline: `6dc80b47fa` (same implementation as `54c146b2b2`).

Aggregate full-AMVF time on 10 ART events falls from 1241.945 to 1218.763 ms:
**1.87% less time (1.019×)**. One-off seeding changes by -0.26%. All seed and
full vertex/track outputs match exactly. Peak process RSS on event 3912208 is
31,384 KiB baseline and 32,256 KiB optimized, an 872 KiB increase.

The track cache now retains the momentum normal-matrix inverse and the weight
matrix after eliminating momentum coordinates. These depend only on the track
covariance, fit dimension and momentum Jacobian. The cache checks covariance
and Jacobian bits and dimension before reuse; shared caches detach on mutation.
The original matrix expressions and their multiplication grouping are retained.
Vertex-position, vertex-covariance and annealing-weight-dependent calculations
still run for every update.

New tests compare cached matrices with direct calculations in 3D and 4D, with
time correlations, repeated calls, dimension changes, covariance changes,
Jacobian-only changes and copied-track independence. All four vertexing suites
pass. Run-3/Run-4 ODD and dense-z stress dumps are byte-identical.

CPU 0, one warmup, five ART samples per event, alternating A/B order. Public
object layouts remain unchanged; only private allocated cache data grows.
These are standalone 2 T replay measurements with the prior ART provenance.
Raw results: [amvf_eleventh_pass_results.json](amvf_eleventh_pass_results.json).
