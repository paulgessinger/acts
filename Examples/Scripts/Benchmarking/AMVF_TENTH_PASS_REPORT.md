# Tenth pass: inline seed refinement helpers (rejected)

Baseline: `54c146b2b2`. Moving updateMaximum and stepSize definitions before
use and marking them inline removes their PLT calls in the compiled object.

Two full 10-event ART blocks show only 0.44% and 0.03% lower AMVF time;
one-off seeding changes by -0.22% and -0.10%. The repeat does not establish a
useful speedup. The candidate was reverted, despite eliminating the calls.
Both blocks have exact outputs and the four targeted suites pass.

CPU 0, one warmup, 5/9 samples per event, alternating A/B order. Data:
[amvf_tenth_pass_results.json](amvf_tenth_pass_results.json).
