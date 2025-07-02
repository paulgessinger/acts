# Blueprint API Guide

The Blueprint system is Acts' Gen3 geometry construction framework that provides a tree-based, declarative approach to building tracking geometries. This guide covers all major Blueprint node types and their usage patterns.

## Core Concepts

### The Blueprint Tree Structure

A Blueprint tree consists of nodes where each node represents a component of the geometry hierarchy:

```
Blueprint (Root)
├── Container Nodes (organize child volumes)
│   ├── Static Volume Nodes (concrete geometry)
│   ├── Layer Nodes (with surfaces)
│   └── Other Container Nodes
├── Material Designator Nodes (material regions)
│   └── Child nodes (inherit material designation)
└── Geometry ID Nodes (identification assignment)
    └── Child nodes (inherit ID configuration)
```

### Construction Phases

All Blueprint construction happens in three phases:

1. **Build Phase**: Volume sizing and bounds calculation
2. **Connect Phase**: Portal creation and inter-volume connections  
3. **Finalize Phase**: TrackingGeometry assembly and navigation setup

## Node Types Reference

### Root Node: Blueprint

The `Blueprint` class is always the root node and entry point for construction.

```cpp
#include "Acts/Geometry/Blueprint.hpp"

using namespace Acts::Experimental;
using namespace Acts::UnitLiterals;

// Create root with world volume envelope
Blueprint::Config config;
config.envelope = Acts::ExtentEnvelope{{
    .r = {10_mm, 15_mm},     // radial envelope [inner, outer]
    .z = {20_mm, 20_mm},     // z envelope [negative, positive]
    .x = {5_mm, 5_mm},       // x envelope [negative, positive] 
    .y = {5_mm, 5_mm}        // y envelope [negative, positive]
}};

auto blueprint = std::make_unique<Blueprint>(config);

// Construct final geometry
Acts::GeometryContext gctx;
auto trackingGeometry = blueprint->construct({}, gctx);
```

### Container Nodes

Container nodes organize child volumes and handle automatic sizing.

#### CylinderContainerBlueprintNode

Arranges child volumes in cylindrical coordinates with automatic sizing.

```cpp
// Add cylinder container with Z-direction stacking
auto& detector = blueprint->addCylinderContainer("Detector", Acts::AxisDirection::AxisZ);

// Configure attachment strategy (how volumes connect)
detector.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::Gap);
detector.setResizeStrategy(Acts::CylinderVolumeStack::ResizeStrategy::Expand);

// R-direction stacking for barrel layers
auto& barrel = detector.addCylinderContainer("Barrel", Acts::AxisDirection::AxisR);
```

**Attachment Strategies:**
- `First`: Attach to first volume
- `Last`: Attach to last volume  
- `Gap`: Insert gap volumes between components
- `Midpoint`: Attach at midpoint

**Resize Strategies:**
- `Expand`: Expand to fit all children
- `Gap`: Insert gaps to fit target size

#### CuboidContainerBlueprintNode

Arranges child volumes in rectangular coordinates.

```cpp
// Create cuboid container
auto& calorimeter = blueprint->addCuboidContainer("ECAL", Acts::AxisDirection::AxisX);

// Add layers in X direction
for (int i = 0; i < 10; ++i) {
    double x = i * 10_mm;
    auto bounds = std::make_shared<Acts::CuboidVolumeBounds>(5_mm, 100_mm, 100_mm);
    
    calorimeter.addStaticVolume(
        Acts::Transform3::Identity() * Acts::Translation3(Acts::Vector3(x, 0, 0)),
        bounds, "ECALLayer" + std::to_string(i));
}
```

### Leaf Nodes

#### StaticBlueprintNode

Represents concrete volumes with fixed geometry.

```cpp
// Method 1: Add with transform and bounds
auto bounds = std::make_shared<Acts::CylinderVolumeBounds>(100_mm, 120_mm, 200_mm);
auto& volume = parent.addStaticVolume(
    Acts::Transform3::Identity(), bounds, "MyVolume");

// Method 2: Add pre-constructed TrackingVolume
auto trackingVol = std::make_unique<Acts::TrackingVolume>(
    Acts::Transform3::Identity(), bounds, "PrebuiltVolume");
auto& volume2 = parent.addStaticVolume(std::move(trackingVol));
```

#### LayerBlueprintNode

Creates layers with surfaces for tracking.

```cpp
// Create layer with surfaces
auto& layer = parent.addLayer("PixelLayer", [&](LayerBlueprintNode& layerNode) {
    // Configure layer properties
    layerNode.setLayerType(LayerBlueprintNode::LayerType::Cylinder)
             .setSurfaces(pixelSurfaces)
             .setTransform(Acts::Transform3::Identity())
             .setEnvelope(Acts::ExtentEnvelope{{.z = {0.5_mm, 0.5_mm}}});
    
    // Optional: Set binning for surface lookup
    layerNode.setBinning(Acts::BinningType::equidistant, Acts::BinningValue::binPhi, 16)
             .setBinning(Acts::BinningType::equidistant, Acts::BinningValue::binZ, 20);
});
```

**Layer Types:**
- `Cylinder`: Cylindrical layer geometry
- `Disc`: Disc-shaped layer geometry  
- `Plane`: Planar layer geometry

### Specialization Nodes

#### MaterialDesignatorBlueprintNode

Designates material regions and inherits to child volumes.

```cpp
// Create material region
auto& materialRegion = parent.addMaterial("SiliconRegion", [&](auto& material) {
    // All children will inherit this material designation
    material.addStaticVolume(transform1, bounds1, "Volume1");
    material.addStaticVolume(transform2, bounds2, "Volume2");
});
```

#### GeometryIdentifierBlueprintNode

Handles geometry identification assignment.

```cpp
// Add geometry identifier node
auto& geoId = parent.withGeometryIdentifier([&](auto& geoIdNode) {
    // Configure ID assignment strategy
    geoIdNode.addStaticVolume(transform, bounds, "IdentifiedVolume");
});
```

## Advanced Examples

### Complex Detector with Multiple Subsystems

```cpp
std::unique_ptr<Acts::TrackingGeometry> buildComplexDetector() {
    using namespace Acts::UnitLiterals;
    
    // Root blueprint
    Blueprint::Config cfg;
    cfg.envelope = Acts::ExtentEnvelope{{.r = {5_mm, 10_mm}, .z = {50_mm, 50_mm}}};
    auto root = std::make_unique<Blueprint>(cfg);
    
    // Inner detector container
    auto& innerDetector = root->addCylinderContainer("InnerDetector", Acts::AxisDirection::AxisZ);
    
    // Pixel detector
    auto& pixel = innerDetector.addCylinderContainer("Pixel", Acts::AxisDirection::AxisZ);
    buildPixelDetector(pixel);
    
    // Strip detector  
    auto& strip = innerDetector.addCylinderContainer("Strip", Acts::AxisDirection::AxisZ);
    buildStripDetector(strip);
    
    // Transition radiation tracker
    auto& trt = innerDetector.addMaterial("TRT_MaterialRegion", [&](auto& material) {
        auto& trtContainer = material.addCylinderContainer("TRT", Acts::AxisDirection::AxisZ);
        buildTRTDetector(trtContainer);
    });
    
    return root->construct({}, Acts::GeometryContext{});
}

void buildPixelDetector(CylinderContainerBlueprintNode& pixel) {
    // Barrel
    auto& barrel = pixel.addCylinderContainer("PixelBarrel", Acts::AxisDirection::AxisR);
    barrel.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::Gap);
    
    for (int layer = 0; layer < 4; ++layer) {
        double r = 50_mm + layer * 40_mm;
        auto& layerNode = barrel.addLayer("B" + std::to_string(layer), [&](auto& ln) {
            ln.setLayerType(LayerBlueprintNode::LayerType::Cylinder)
              .setSurfaces(createPixelBarrelSurfaces(r))
              .setEnvelope(Acts::ExtentEnvelope{{.r = {2_mm, 2_mm}, .z = {5_mm, 5_mm}}});
        });
    }
    
    // Positive endcap
    auto& posEC = pixel.addCylinderContainer("PixelPosEndcap", Acts::AxisDirection::AxisZ);
    for (int disc = 0; disc < 3; ++disc) {
        double z = 400_mm + disc * 150_mm;
        auto& layerNode = posEC.addLayer("EC+_" + std::to_string(disc), [&](auto& ln) {
            ln.setLayerType(LayerBlueprintNode::LayerType::Disc)
              .setSurfaces(createPixelEndcapSurfaces(z))
              .setTransform(Acts::Transform3::Identity() * 
                           Acts::Translation3(Acts::Vector3(0, 0, z)));
        });
    }
    
    // Negative endcap  
    auto& negEC = pixel.addCylinderContainer("PixelNegEndcap", Acts::AxisDirection::AxisZ);
    for (int disc = 0; disc < 3; ++disc) {
        double z = -400_mm - disc * 150_mm;
        auto& layerNode = negEC.addLayer("EC-_" + std::to_string(disc), [&](auto& ln) {
            ln.setLayerType(LayerBlueprintNode::LayerType::Disc)
              .setSurfaces(createPixelEndcapSurfaces(z))
              .setTransform(Acts::Transform3::Identity() * 
                           Acts::Translation3(Acts::Vector3(0, 0, z)));
        });
    }
}
```

### Working with Python Blueprints

The Blueprint system is fully exposed to Python for interactive development:

```python
import acts

# Create root blueprint
root = acts.Blueprint(envelope=acts.ExtentEnvelope(r=[10*acts.mm, 10*acts.mm]))

# Use context managers for clean hierarchy
with root.CylinderContainer("Detector", acts.AxisDirection.AxisZ) as detector:
    detector.attachmentStrategy = acts.CylinderVolumeStack.AttachmentStrategy.Gap
    
    with detector.CylinderContainer("Barrel", acts.AxisDirection.AxisR) as barrel:
        for i in range(4):
            r = 50 + i * 30  # mm
            bounds = acts.CylinderVolumeBounds(r*acts.mm, (r+5)*acts.mm, 300*acts.mm)
            barrel.addStaticVolume(acts.Transform3.Identity(), bounds, f"Layer{i}")
    
    with detector.CylinderContainer("Endcap", acts.AxisDirection.AxisZ) as endcap:
        for i in range(3):
            z = 350 + i * 100  # mm
            bounds = acts.CylinderVolumeBounds(100*acts.mm, 200*acts.mm, 25*acts.mm)
            transform = acts.Transform3.Identity() * acts.Translation3(acts.Vector3(0, 0, z*acts.mm))
            endcap.addStaticVolume(transform, bounds, f"Disc{i}")

# Construct and visualize
gctx = acts.GeometryContext()
geometry = root.construct(options=acts.BlueprintNode.Options(), gctx=gctx)

# Export tree structure for debugging
with open("detector_tree.dot", "w") as f:
    root.graphviz(f)
```

## Best Practices

### 1. Use Appropriate Container Types
- **CylinderContainer**: For detectors with cylindrical symmetry (barrels, endcaps)
- **CuboidContainer**: For rectangular calorimeters, muon chambers

### 2. Configure Attachment Strategies
```cpp
// For detector layers that should touch
container.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::First);

// For detector components that need spacing
container.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::Gap);
```

### 3. Use Material and GeometryID Nodes Appropriately
```cpp
// Group related volumes under material designators
auto& siliconRegion = root->addMaterial("Silicon", [&](auto& material) {
    material.addCylinderContainer("PixelDetector", Acts::AxisDirection::AxisZ);
    material.addCylinderContainer("StripDetector", Acts::AxisDirection::AxisZ);
});
```

### 4. Leverage Callbacks for Clean Code
```cpp
// Use callbacks to encapsulate configuration
auto& detector = root->addCylinderContainer("Detector", Acts::AxisDirection::AxisZ, 
    [&](auto& container) {
        container.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::Gap);
        // Add children...
    });
```

### 5. Debug with Tree Visualization
```cpp
// Export tree structure for debugging
std::ofstream dotFile("geometry_tree.dot");
blueprint->graphviz(dotFile);
// Use graphviz to render: dot -Tpng geometry_tree.dot -o tree.png
```

The Blueprint system provides a powerful, intuitive way to construct complex detector geometries while automatically handling the intricate details of portal management, volume sizing, and TrackingGeometry assembly.