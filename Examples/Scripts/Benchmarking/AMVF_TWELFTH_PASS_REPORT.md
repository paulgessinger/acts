# Twelfth pass: derivative batching experiment (rejected)

Baseline: `899657fdeb`. Buffer eight contributing tracks, calculate derivatives
in a fixed-size loop, then add each channel in original track order. Retain
scalar libm exp and a scalar fallback below 32 tracks.

The compiled object contains packed arithmetic, and all ART output comparisons
and four vertexing suites pass. Nevertheless, aggregate full-AMVF time rises
from 1215.358 to 1292.237 ms (**6.33% slower**) and standalone seeding rises
12.19%. Buffering and marshaling do not pay off for this data flow. The change
was reverted; no additional repeat was needed for this clear regression.

CPU 0, one warmup and five samples/event, alternating A/B order. Data:
[amvf_twelfth_pass_results.json](amvf_twelfth_pass_results.json).
