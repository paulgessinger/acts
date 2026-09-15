# Thirteenth pass: interval-bin granularity (retain 64)

Baseline code: `ca24cdae86` / `899657fdeb`, with 64 bins.

| Candidate | Full-AMVF time change | Direct seeder time change |
|---|---:|---:|
| 32 bins | +1.76% | +6.21% |
| 128 bins | +1.28% | -3.45% |
| 256 bins | +7.33% | -3.78% |

The larger indices help one-off seeding but lose in repeated full-AMVF searches.
None improves the main workload, so all candidates were reverted and 64 bins
retained. Each variant passed the exhaustive density oracle and exact seed and
full vertex/track comparisons on the 10 ART events.

One warmup and five repetitions/event, CPU 0, alternating A/B order; each
candidate was compared with the same 64-bin baseline library. Raw data:
[amvf_thirteenth_pass_results.json](amvf_thirteenth_pass_results.json).
