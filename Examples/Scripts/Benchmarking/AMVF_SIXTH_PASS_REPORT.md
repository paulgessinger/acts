# Sixth pass: cache insertion hint experiment (rejected)

Baseline: `6e6eb96142`. Replace map find + emplace with lower_bound +
emplace_hint to avoid a second tree search on a cache miss.

Two complete ART A/B blocks (5 then 9 repetitions/event, one warmup, CPU 0)
passed exact seed and full-output checks. The first block had several noisy
measurements and showed a 3.4% full-AMVF regression. The repeat showed only
0.42% lower full-AMVF time, while direct seeding regressed 2.1% despite caching
being disabled for that path. The small full-AMVF gain does not justify keeping
the change and cold-path regression. The source change was reverted.

Four targeted vertexing suites passed. Results are retained in
[amvf_sixth_pass_results.json](amvf_sixth_pass_results.json).
