# Seventh pass: pool Gaussian cache nodes

Baseline code is unchanged between `6e6eb96142` and `5e5955b2c3`.

A per-EvaluationCache unsynchronized pool resource now supplies ordered-map
nodes. Range invalidation returns nodes to the pool for subsequent queries.
The resource belongs to one state, and copied states construct their own pool.
Map ordering, query equality, strict interval invalidation and capacity limits
are unchanged. Pool storage follows peak cache demand during that state lifetime.

Two complete 10-event ART blocks show **1.17% and 1.22% lower full-AMVF time**.
One-off seeding changes by +0.34% and +0.54%. Peak process RSS on event 3912208
(10 measured fits) is 22,016 KiB baseline and 22,180 KiB optimized. This is one
process-level observation, not a general allocation bound.

Both blocks have exact seed and byte-identical vertex/track outputs. Four
vertexing suites pass, including copied-density-state cache tests. Run-3 and
Run-4 ODD/dense stress replays also have byte-identical outputs.

Timings use CPU 0, one warmup and 5/9 measured samples per ART event, alternating
A/B order. Public state layouts are unchanged. These remain standalone 2 T
replays with the prior ART input/beam metadata.

Raw data: [amvf_seventh_pass_results.json](amvf_seventh_pass_results.json).
