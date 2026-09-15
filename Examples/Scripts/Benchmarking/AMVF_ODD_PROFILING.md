# AMVF ODD baseline and profiling harness

`ActsBenchmarkAdaptiveMultiVertexFinder` isolates the
`AdaptiveMultiVertexFinder` and its `AdaptiveMultiVertexFitter` from the rest
of reconstruction. It replays fitted track parameters and their diagonal
covariance estimates exported from an ODD `tracksummary_ambi.root` file.
Simulation, digitisation, CKF, ambiguity resolution, Athena track selection,
input parsing, and setup are outside of the timed region.

The default `--athena-era run4` profile mirrors the effective ordinary offline
primary-vertex configuration on Athena `main`: Gaussian density seeding, the
Athena component defaults, the Run-4 `minWeight=0.02` and
`maxIterations=200` overrides, the ITk `tracksMaxZinterval=0.5 mm`, and an
enabled beam-spot constraint. `--athena-era run3` instead selects the ID
`tracksMaxZinterval=3 mm` and the component defaults `minWeight=0.0001` and
`maxIterations=100`. Both profiles use `EigenStepper`, as does Athena.

The default zero-centred beam spot with widths `(0.15, 0.15, 53) mm` is the
fallback in Athena's `BeamSpotCondAlg`. Production obtains these quantities
from conditions data, so use the `--beamspot-sigma-{x,y,z}` options when
replaying an event with known IOV values.

## Representative input

Generate a deterministic ODD particle-gun event with 200 primary vertices:

```console
python Examples/Scripts/Python/full_chain_odd.py \
  --events 1 --jobs 1 \
  --gun-multiplicity 200 --gun-particles 4 \
  --output-root --no-output-csv \
  --output odd-mu200
```

This is a controlled Run-4-complexity proxy, not a claim that four muons per
vertex reproduce the full track spectrum of minimum-bias pile-up. For final
optimization numbers, replay an ODD `ttbar + mu=200` sample produced with
`--ttbar --ttbar-pu 200`; the harness is unchanged.

Export the six bound parameters and their uncertainties to a flat CSV once.
The expected header is:

```text
loc0,loc1,phi,theta,qop,time,sigma_loc0,sigma_loc1,sigma_phi,sigma_theta,sigma_qop,sigma_time
```

The included ROOT macro performs this export without putting ROOT I/O in the
measured executable:

```console
root -l -b -q \
  'Examples/Scripts/Benchmarking/export_amvf_tracks.C("odd-mu200/tracksummary_ambi.root","odd-mu200-tracks.csv")'
```

## Timing baseline

Configure a symbol-bearing build from the benchmark worktree and run:

```console
cmake -S . -B build-amvf -GNinja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DACTS_BUILD_BENCHMARKS=ON
cmake --build build-amvf --target ActsBenchmarkAdaptiveMultiVertexFinder

taskset -c 0 build-amvf/bin/ActsBenchmarkAdaptiveMultiVertexFinder \
  --input odd-mu200-tracks.csv --athena-era run4 \
  --warmup 1 --repetitions 7 \
  | tee amvf-odd-baseline.json
```

The stable comparison number is `median_ms`. Keep the input CSV, compiler,
build type, CPU affinity, and benchmark configuration fixed when comparing
changes.

## CPU profile

Use a `RelWithDebInfo` build and profile several repetitions to collect enough
samples:

```console
perf record -F 199 --call-graph dwarf \
  -o amvf-perf.data -- \
  taskset -c 0 build-amvf/bin/ActsBenchmarkAdaptiveMultiVertexFinder \
  --input odd-mu200-tracks.csv --athena-era run4 \
  --warmup 1 --repetitions 5

perf report --stdio --no-children \
  -i amvf-perf.data > amvf-perf-flat.txt

perf report --stdio --children \
  -i amvf-perf.data > amvf-perf-callgraph.txt
```

The executable performs input parsing and object construction before the
repeated fit, so these one-off costs appear in `perf` but not in `median_ms`.
Use the flat report's ACTS vertexing symbols to rank optimization targets.

## Comparing reconstruction results

Use `--output-vertices vertices.txt` to save the final repetition outside the
measured region. The output includes vertex positions, covariances, fit quality,
input-track indices, weights, compatibility, and fitted track parameters and
covariances at 17-digit precision. Compare the files with `cmp` for identical
inputs and configuration. This detects changes that a vertex-count comparison
would miss.

See [the optimization report](AMVF_ODD_OPTIMIZATION_REPORT.md) for the measured
support-interval lookup optimization and its validation.
