# AMVF fifth pass: cache input global z positions

Baseline: `5f72399c3b`.

On the 10 Run-4 ART events, aggregate full-AMVF time (sum of per-event medians)
falls from 1572.585 to 1375.893 ms: **12.51% less time, 1.143× speedup**.
All seed results compare exactly and all full vertex/track dumps are byte-identical.

The finder calculates input global z positions once per find invocation using
the configured parameter extractor, reference surfaces and geometry context.
Compatible-track searches and nearest-track recovery reuse these coordinates.
The cache is local to the call and uses input identity; it cannot retain values
from a previous event/context. Tracks retain their original traversal order and
the strict z-window predicate is unchanged. Input parameters are treated as
immutable during find, as elsewhere in the vertexing algorithms. Candidates
passing the window still use the original impact-significance calculation.

Measurement uses CPU 0, one warmup and five samples per event, alternating A/B
process order. Each variant uses the same executable with its matching library;
public state layouts are unchanged. Four targeted suites pass, including the
user-defined track and 4D fitter tests.

Raw data: [amvf_fifth_pass_results.json](amvf_fifth_pass_results.json).
The same ART provenance and constant-2-T standalone propagation limitation
apply; see [the sample report](AthenaExport/RESULTS.md).

Run-3/Run-4 ODD and dense-z stress comparisons also have byte-identical outputs.
A transformed-reference ART case (rotation 0.03 rad around y and translation
(0.01,-0.02,2) mm) passes exact A/B checks in both eras. This exercises global
position conversion beyond the original origin perigee surfaces.
