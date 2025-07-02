# Geometry Construction Approaches: Gen1 vs Gen3

Acts provides two distinct approaches for constructing tracking geometries. The original "Gen1" approach uses traditional builders and creators, while the newer "Gen3" approach introduces a tree-based Blueprint system that provides a more declarative and hierarchical way to construct geometries.

## Overview

### Gen1 (Original) Approach
The Gen1 approach uses a collection of specialized builder classes that must be manually orchestrated to create geometry:
- `LayerCreator` - Creates individual layers from surfaces
- `CylinderVolumeHelper` - Helps construct cylindrical volumes
- `LayerArrayCreator` - Creates arrays of layers
- `TrackingVolumeArrayCreator` - Creates arrays of tracking volumes
- `SurfaceArrayCreator` - Creates surface arrays for efficient lookup

### Gen3 (Blueprint) Approach  
The Gen3 approach uses a tree-based Blueprint system where:
- `Blueprint` is the root node that coordinates construction
- `BlueprintNode` is the base class for all tree nodes
- Specialized nodes handle different construction aspects
- Construction happens in three coordinated phases: Build, Connect, Finalize

## Comparison

| Aspect | Gen1 | Gen3 |
|--------|------|------|
| **Paradigm** | Imperative, manual assembly | Declarative, tree-based |
| **Code Style** | Procedural with builders | Hierarchical with nodes |
| **Complexity** | High - requires detailed knowledge | Lower - intuitive structure |
| **Maintainability** | Difficult - scattered logic | Better - organized hierarchy |
| **Debugging** | Hard to trace construction flow | Clear tree structure |
| **Flexibility** | Limited by builder interfaces | Highly composable nodes |

## Gen1 Example: Creating a Pixel Detector

```cpp
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"

// Gen1 approach - manual construction
void buildPixelDetectorGen1(const GeometryContext& gctx) {
    using namespace Acts::UnitLiterals;
    
    // 1. Create individual components
    Acts::LayerCreator::Config layerConfig;
    layerConfig.surfaceArrayCreator = std::make_shared<Acts::SurfaceArrayCreator>();
    Acts::LayerCreator layerCreator(layerConfig, logger);
    
    Acts::CylinderVolumeHelper::Config volConfig;
    volConfig.layerArrayCreator = std::make_shared<Acts::LayerArrayCreator>();
    volConfig.trackingVolumeArrayCreator = std::make_shared<Acts::TrackingVolumeArrayCreator>();
    Acts::CylinderVolumeHelper volumeHelper(volConfig, logger);
    
    // 2. Manually create surfaces (not shown - complex process)
    std::vector<std::shared_ptr<const Acts::Surface>> barrelSurfaces = createBarrelSurfaces();
    std::vector<std::shared_ptr<const Acts::Surface>> endcapSurfaces = createEndcapSurfaces();
    
    // 3. Create layers from surfaces
    auto barrelLayer = layerCreator.cylinderLayer(gctx, barrelSurfaces, 16, 20);
    auto endcapLayer = layerCreator.discLayer(gctx, endcapSurfaces, 10, 16);
    
    // 4. Assemble layers into volumes
    Acts::LayerVector barrelLayers = {barrelLayer};
    Acts::LayerVector endcapLayers = {endcapLayer};
    
    auto barrelBounds = std::make_shared<Acts::CylinderVolumeBounds>(50_mm, 200_mm, 300_mm);
    auto endcapBounds = std::make_shared<Acts::CylinderVolumeBounds>(100_mm, 200_mm, 50_mm);
    
    auto barrelVolume = volumeHelper.createTrackingVolume(
        gctx, barrelLayers, barrelBounds, Acts::Transform3::Identity(), "PixelBarrel");
    
    auto endcapVolume = volumeHelper.createTrackingVolume(
        gctx, endcapLayers, endcapBounds, 
        Acts::Transform3::Identity() * Acts::Translation3(Acts::Vector3(0, 0, 400_mm)),
        "PixelEndcap");
    
    // 5. Create container volume
    std::vector<Acts::TrackingVolumePtr> volumes = {barrelVolume, endcapVolume};
    auto containerBounds = std::make_shared<Acts::CylinderVolumeBounds>(50_mm, 250_mm, 500_mm);
    
    auto containerVolume = volumeHelper.createContainerTrackingVolume(
        gctx, volumes, containerBounds, Acts::Transform3::Identity(), "PixelDetector");
}
```

## Gen3 Example: Creating the Same Pixel Detector

```cpp
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"

// Gen3 approach - declarative tree construction
std::unique_ptr<Acts::TrackingGeometry> buildPixelDetectorGen3(const GeometryContext& gctx) {
    using namespace Acts::UnitLiterals;
    using namespace Acts::Experimental;
    
    // 1. Create root blueprint with envelope
    Blueprint::Config cfg;
    cfg.envelope = Acts::ExtentEnvelope{{
        .r = {5_mm, 10_mm},    // radial envelope
        .z = {10_mm, 10_mm}    // z envelope  
    }};
    auto root = std::make_unique<Blueprint>(cfg);
    
    // 2. Add pixel detector container
    auto& pixel = root->addCylinderContainer("Pixel", Acts::AxisDirection::AxisZ);
    
    // 3. Add barrel section
    auto& barrel = pixel.addCylinderContainer("PixelBarrel", Acts::AxisDirection::AxisR);
    barrel.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::Gap);
    
    // Add barrel layers with automatic sizing
    for (int i = 0; i < 4; ++i) {
        double radius = 50_mm + i * 30_mm;
        auto bounds = std::make_shared<Acts::CylinderVolumeBounds>(
            radius, radius + 5_mm, 300_mm);
        
        barrel.addStaticVolume(
            Acts::Transform3::Identity(), bounds, 
            "PixelBarrelLayer" + std::to_string(i));
    }
    
    // 4. Add positive endcap
    auto& posEndcap = pixel.addCylinderContainer("PixelPosEndcap", Acts::AxisDirection::AxisZ);
    posEndcap.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::Gap);
    
    for (int i = 0; i < 3; ++i) {
        double z = 350_mm + i * 100_mm;
        auto bounds = std::make_shared<Acts::CylinderVolumeBounds>(100_mm, 180_mm, 25_mm);
        
        posEndcap.addStaticVolume(
            Acts::Transform3::Identity() * Acts::Translation3(Acts::Vector3(0, 0, z)),
            bounds, "PixelPosEndcapDisk" + std::to_string(i));
    }
    
    // 5. Add negative endcap (with material designation)
    auto& material = pixel.addMaterial("PixelMaterialRegion");
    auto& negEndcap = material.addCylinderContainer("PixelNegEndcap", Acts::AxisDirection::AxisZ);
    
    for (int i = 0; i < 3; ++i) {
        double z = -350_mm - i * 100_mm;
        auto bounds = std::make_shared<Acts::CylinderVolumeBounds>(100_mm, 180_mm, 25_mm);
        
        negEndcap.addStaticVolume(
            Acts::Transform3::Identity() * Acts::Translation3(Acts::Vector3(0, 0, z)),
            bounds, "PixelNegEndcapDisk" + std::to_string(i));
    }
    
    // 6. Construct the geometry (automatic 3-phase construction)
    Acts::BlueprintOptions options;
    return root->construct(options, gctx);
}
```

## Key Advantages of Gen3

### 1. Declarative Structure
The tree structure clearly shows the hierarchy:
```
Pixel (Container)
├── PixelBarrel (Container)
│   ├── PixelBarrelLayer0 (Static Volume)
│   ├── PixelBarrelLayer1 (Static Volume)
│   └── ...
├── PixelPosEndcap (Container)
│   ├── PixelPosEndcapDisk0 (Static Volume)
│   └── ...
└── PixelMaterialRegion (Material Designator)
    └── PixelNegEndcap (Container)
        └── ...
```

### 2. Automatic Portal Management
Gen3 automatically handles portal creation and connection between volumes, eliminating manual portal management that was error-prone in Gen1.

### 3. Intelligent Sizing
Container nodes automatically size themselves based on their children, with configurable attachment and resize strategies.

### 4. Composable Design
Nodes can be easily combined and reused. Material designation, geometry identification, and other aspects are handled by specialized nodes.

## Tree Construction Phases

The Gen3 system uses a three-phase construction process:

### Phase 1: Build
- Each node builds its volume representation
- Size calculation and propagation up the tree
- Container nodes size themselves based on children

### Phase 2: Connect  
- Portal creation at volume boundaries
- Portal connection between adjacent volumes
- Portal shell creation for complex geometries

### Phase 3: Finalize
- Portal registration with volumes
- Navigation policy creation
- Geometry identifier assignment
- Final assembly into TrackingGeometry

## Migration Guidelines

When migrating from Gen1 to Gen3:

1. **Identify the hierarchy** - Map your existing volumes to a tree structure
2. **Choose container types** - Use `CylinderContainer` for cylindrical arrangements, `CuboidContainer` for rectangular
3. **Replace builders with nodes** - Convert `LayerCreator` calls to `addLayer()` nodes
4. **Use automatic sizing** - Let containers handle size calculation instead of manual bounds computation
5. **Leverage specialization nodes** - Use `MaterialDesignatorBlueprintNode` for material regions, `LayerBlueprintNode` for surfaces

The Gen3 approach significantly reduces boilerplate code while providing better maintainability and clearer structure for complex detector geometries.