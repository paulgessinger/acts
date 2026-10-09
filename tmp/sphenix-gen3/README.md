# sPHENIX: Gen3 geometry from the TGeo file (discussion only, not for merging)

Builds a Gen3 (blueprint) tracking geometry from `sPHENIXActsGeom.root` using
the TGeo `BlueprintBuilder` (`ActsPlugins/Root/BlueprintBuilder.hpp`), trying to
resemble what the Gen1 `TGeoDetector` + `tgeo-sphenix-mms-actsv45.0.0.json`
produce. There are no Python bindings for the TGeo blueprint builder, so this is
a small standalone C++ program.

## Findings that shaped it

- **MVTX is off the beam axis** by about (6.3, -1.1, -0.7) mm
  (`log_MVTX_Wrapper` placement). Seen from the beam line its three layers
  overlap in r, so they can't be separate coaxial layer volumes. Gen1 also builds
  a single MVTX layer (`geo-tgeo-layer-r-split: 0`).
- **INTT has 4 ladder layers** (`ladder_<layer>_...`), but each barrel's two
  staggered sub-layers overlap radially, so at most 2 layer volumes are possible.
  Gen1's layers 4 and 8 in volume 12 (the cut-off approach surfaces near
  eta = 0.7) come from a bug in `ProtoLayerHelper`'s greedy clustering: a sensor
  that bridges two clusters only joins the first one, leaving small stray
  "layers" with short approach surfaces.
- **TPC measurement volumes are ~5.6 mm thick gas slabs.** Surface thickness
  enters the layer extent, so neighbouring TPC layers overlap (Gen1's approach
  surfaces overlap too). Gen3 volumes can't, so TPC surface thickness is set to 0.

## Layout (`sphenix_gen3.cpp`)

| Subsystem  | Volumes          | Navigation    | Material surfaces                                     |
|------------|------------------|---------------|-------------------------------------------------------|
| MVTX       | 1                | TryAll        | 3 passive carrier cylinders, centred on the MVTX axis |
| INTT       | 1                | TryAll        | 4 passive carriers, one per ladder layer              |
| TPC        | 48 layer volumes | surface array | layer portals (inner + outer)                         |
| Micromegas | 1                | TryAll        | 2 passive carriers (phi tiles, z tiles)               |

Carriers are cylinder surfaces inside the volume with `ProtoGridSurfaceMaterial`
(36 x 50 rphi-z bins, deferred axes). Their r / z come from the sensors, so they
are independent of how the volumes are split, unlike Gen1 approach surfaces.

`--zgaps` wraps each subsystem in a z-stack with gap volumes, reproducing Gen1's
`fGap | Barrel | sGap` structure (MVTX +-137 mm, INTT +-232 mm, TPC +-1027 mm).
Without it, the R-stack expands every volume to the full world length.

A 5000-track straight-line navigation check (|eta| < 1.2) runs at the end:
0 failures, about 3 MVTX / 2 INTT / 47.5 TPC sensitive hits and 3 / 4 MVTX / INTT
carrier crossings per track.

Open: material *mapping* onto the internal carriers hasn't been tested yet.

## Plots

- `sphenix_gen3_rz.png`: r-z, full detector and inner silicon, both z modes
- `sphenix_gen3_xy.png`: x-y, full detector, MVTX + INTT, MVTX zoom with axis
- `gen1_vs_gen3_rz_layered.png`: Gen1 approach surfaces vs. an earlier,
  layer-volume-only Gen3 variant (`sphenix_gen3_layered.cpp`)

## Running

```sh
ACTS_BUILD=/path/to/acts/build ./build.sh            # needs ACTS_BUILD_PLUGIN_ROOT=ON
./sphenix_gen3 sPHENIXActsGeom.root out/gen3 [--zgaps] [loglevel]
python plot_layout.py out                            # expects out/gen3_* and out/gen3_zgaps_*
python gen1_csv.py sPHENIXActsGeom.root tgeo-sphenix-mms-actsv45.0.0.json out/gen1_approach.csv
```

The program writes `<prefix>_volumes.csv`, `_sensors.csv`,
`_material_surfaces.csv`, `_hits.csv` and `_blueprint.dot` (graphviz).
