(geometry_impl)=

# Geometry module

The ACTS geometry model is strongly based on the ATLAS Tracking geometry. Its
core is built on a surface-based description that make up all geometry objects
of higher complexity. This design has been chosen as the surface objects can be
used together with the track propagation module and thus all geometry objects
become natively integrated into the tracking software.

```{note}
There is an ongoing rewrite of the geometry and navigation modules where
logical layers will be modelled as volumes, see [](layerless_geometry).

The Acts geometry system now provides two construction approaches:
- **Gen1 (Original)**: Traditional builders and manual assembly
- **Gen3 (Blueprint)**: Tree-based declarative construction with automatic portal management

For new projects, use the Gen3 Blueprint system. See the migration guide for transitioning existing code.

```

:::{toctree}
:maxdepth: 1
concepts
construction_approaches
blueprint_api_guide
blueprint_diagrams
gen1_to_gen3_migration
geometry_id
material
surfaces
legacy/legacy
construction
layerless/layerless
:::
