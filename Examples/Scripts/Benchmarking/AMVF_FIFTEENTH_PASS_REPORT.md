# Fifteenth pass: contiguous scans for full bins (rejected)

Baseline: `a60a7042f9`, with 64 bins and the 8 MiB exponential table.

The candidate returned the exhaustive contiguous scan whenever an interval bin
contained every density entry. This avoids indirect access for that bin while
preserving the original bounds, arithmetic and accumulation order.

Full-AMVF time on the ten tuning ART events changed by **-0.06%** and direct
seeding by **+0.04%**. This is no useful gain. The dense-z Run-3 stress case
became 4.45% slower, and the remaining stress checks did not establish a reason
to keep the extra branch. The candidate was reverted.

All four vertexing suites passed. Exact seed and full vertex/track comparisons
passed for all ten ART events and the Run-3/Run-4 ODD, dense-z and transformed
perigee fixtures. One warmup and five repetitions, CPU 0; ART A/B order alternated
by event. Reconstruction had finished before timing began.

Raw measurements: [amvf_fifteenth_pass_results.json](amvf_fifteenth_pass_results.json).
