# Athena ART track export and standalone AMVF replay

This small Athena package exports the track population selected by the production
vertex-track selector, without changing the reconstruction algorithms. It writes
one full-covariance ACTS `TrackParameterData` CSV and one metadata file per event.
Upstream reconstruction is run once; subsequent timing needs only these files.

## Build the exporter

In a clean Bash shell, from the ACTS worktree:

```sh
source /cvmfs/sft.cern.ch/lcg/releases/gcc/15.2.0/x86_64-el9/setup.sh
source /cvmfs/atlas-nightlies.cern.ch/repo/sw/main_Athena_x86_64-el9-gcc15-opt/2026-09-14T2100/Athena/25.0.73/InstallArea/x86_64-el9-gcc15-opt/setup.sh
/usr/bin/cmake -S Examples/Scripts/Benchmarking/AthenaExport/Projects/WorkDir \
  -B /tmp/amvf-athena-export -GNinja
/usr/bin/cmake --build /tmp/amvf-athena-export -j4
```

Use the native CMake executable: a Python-installed CMake wrapper can conflict
with Athena's Python environment. The exporter is a separate WorkDir package;
it does not modify the Athena checkout or rebuild ACTS inside Athena.

## Produce the sample

In a fresh shell:

```sh
Examples/Scripts/Benchmarking/AthenaExport/run_export.sh \
  /absolute/path/to/new-output-directory 10
```

The script uses the production-flags leg of the Run-4 PU200 ART test, with
`OnlyTrackingRecoPreInclude` and `actsProductionFlags`. It resolves input and
conditions from the pinned release's `TestDefaults`, saves those exact values
in `export-provenance.json`, and records the evaluated configuration flags and
configuration pickle. It uses normal Frontier conditions access and the CVMFS
conditions-file catalog; network access to Frontier is required.

`AMVF_ATHENA_RELEASE` and `AMVF_EXPORT_BUILD` can select a different release/build
pair. The exporter must be built against that release. The reference-data path
can be redirected through the usual `ATLAS_REFERENCE_DATA` environment variable.
`AMVF_TRACK_COLLECTION` defaults to `InDetTrackParticles` for the pure ACTS chain.
Use a fresh output directory so exports from different runs cannot be mixed.

## Data contract

- `eventNNNNNNNNN-tracks.csv`: original track index, five perigee parameters,
  all 25 spatial covariance elements. This is the existing ACTS CSV schema.
- `eventNNNNNNNNN-metadata.txt`: schema version, run/event/lumiblock, average and
  actual mu, input/selected/missing-covariance counts, beam-constraint flag,
  beam position, full 3x3 beam covariance, and 4x4 perigee transform. Matrices
  are row-major; floating-point numbers are written at 17-digit precision.
- Units are mm and GeV: q/p is multiplied by 1000, and its covariance row and
  column are each multiplied by 1000. The no-time replay initializes t=0,
  Ctt=1 and all time correlations to zero, matching the Athena vertex wrapper.

The exporter calls `VtxInDetTrackSelectionCfg` with the same flags and beam
reference as the vertex tool. It then reads `perigeeParameters()`, preserving
the conversion's full covariance. It preserves input order and checks that
selected tracks share a reference surface, as assumed by Athena's wrapper.
Beam-spot information is read from the actual `BeamSpotData` conditions object.

## Verify and replay

Build `ActsBenchmarkAdaptiveMultiVertexFinder` as described in the parent
profiling document. The benchmark reads the standard CSV schema directly to
keep this executable independent of the examples sequencer. The file is also
compatible with `CsvTrackParameterReader`; its perigee configuration must be
set from the metadata. The benchmark additionally accepts the full transform
and beam constraint through `--metadata`.

```sh
python3 Examples/Scripts/Benchmarking/AthenaExport/verify_export.py \
  /tmp/acts-amvf-optimized-build/bin/ActsBenchmarkAdaptiveMultiVertexFinder \
  /absolute/path/to/new-output-directory/tracks

python3 Examples/Scripts/Benchmarking/AthenaExport/replay_exports.py \
  /tmp/acts-amvf-optimized-build/bin/ActsBenchmarkAdaptiveMultiVertexFinder \
  /absolute/path/to/new-output-directory/tracks \
  --baseline-library /tmp/acts-amvf-build-main/lib64 \
  --optimized-library /tmp/acts-amvf-optimized-build/lib64 \
  --repetitions 5 --output /absolute/path/to/replay-results.json
```

The verifier compares every parsed parameter and covariance element against the
CSV. The replay driver alternates A/B order across events, pins the executable
to one CPU, checks exact seed position/width equality, and compares vertex/track
output dumps byte-for-byte. Input parsing and output writing are outside timing.
Library substitution requires compatible public class layouts. When layouts
change, build each revision's executable against its own headers and pass
`--baseline-binary /path/to/baseline/ActsBenchmarkAdaptiveMultiVertexFinder`.
The driver then runs each executable with its corresponding library. This is
required for the third pass, which extends the density/finder state.

The seeder-only measurement uses the production track covariance directly and
requires no field or detector geometry. Full-AMVF replay retains the standalone
2 T constant field and simple propagation setup. It is an A/B comparison on
realistic inputs, not an exact replay of Athena's field and navigator.
