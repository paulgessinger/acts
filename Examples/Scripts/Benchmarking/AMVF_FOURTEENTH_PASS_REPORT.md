# Fourteenth pass: expand the scalar Gaussian table to 16 MiB

The retained cache limit is **16 MiB of scalar entries, at most 2048 rows**, up
from 8 MiB / 1024 rows. The larger table covers more initial trial positions
without changing which tracks contribute or the order of arithmetic.

## Measurements

| Comparison | Events | Full-AMVF time change | Direct seeder change |
|---|---:|---:|---:|
| 8 → 16 MiB, retained change | 10 tuning | -2.43% | -0.01% |
| 8 → 32 MiB, exploratory | 10 tuning | -2.78% | +0.03% |
| 8 → 32 MiB, exploratory | 90 additional | -3.53% | -0.49% |
| 16 → 32 MiB | All 100 | -0.22% | -0.38% |

The retained 16 MiB variant reduces the sum of ten per-event AMVF medians from
1216.040 to 1186.488 ms. Doubling again to 32 MiB gives no compelling additional
gain across all 100 events; its 0.22% difference is smaller than the corresponding
cold-seeder variation, whose scalar table is disabled. Keep the smaller limit.
These are separately measured blocks; do not multiply their percentages.

The 32 MiB prototype adds 6532 KiB peak process RSS on event 3912208 and
13224 KiB on the sample's highest-track-count event relative to the 8 MiB version.
The scalar limit excludes row/container overhead and the rest of the event state.
Final 16 MiB memory measurements are recorded in the final validation report.

## Validation and baseline

All four targeted vertexing suites pass. Every seed and full vertex/track dump
matches exactly in all four ART comparison blocks. Run-3/Run-4 ODD, dense-z and
transformed-perigee comparisons also pass for both cache sizes.

The initial 32 MiB comparison uses the kept implementation from `899657fdeb`.
The final 8/16 MiB and 16/32 MiB comparisons both include `a60a7042f9`, which
defers unused position extraction. Each pair differs only in the cache bound.
The rejected full-bin shortcut is absent.

CPU 0, one warmup, alternating A/B order by event; five measurements per event
for the ten-event blocks and three for the larger comparisons. Stress checks
use five repetitions for 32 MiB and three for 16 MiB. Athena reconstruction and
all compilation finished before timing. The existing constant-2-T standalone
replay scope applies.

Raw results: [amvf_fourteenth_pass_results.json](amvf_fourteenth_pass_results.json).
