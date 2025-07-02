# Migration Guide: Gen1 to Gen3 Geometry Construction

This guide provides step-by-step instructions for migrating existing Gen1 geometry construction code to the new Gen3 Blueprint system.

## Quick Migration Checklist

- [ ] Replace `LayerCreator` with `LayerBlueprintNode`
- [ ] Replace `CylinderVolumeHelper` with `CylinderContainerBlueprintNode`
- [ ] Replace `TrackingVolumeArrayCreator` with container nodes
- [ ] Replace manual volume bounds calculation with automatic sizing
- [ ] Replace manual portal management with automatic portal system
- [ ] Update includes to use Blueprint headers
- [ ] Test with new three-phase construction system

## Component Mapping

### Builder Class Replacements

| Gen1 Component | Gen3 Replacement | Purpose |
|---------------|------------------|---------|
| `LayerCreator` | `LayerBlueprintNode` | Creating layers with surfaces |
| `CylinderVolumeHelper` | `CylinderContainerBlueprintNode` | Cylindrical volume arrangement |
| `LayerArrayCreator` | Container node automatic sizing | Layer organization |
| `TrackingVolumeArrayCreator` | Container node automatic sizing | Volume organization |
| `SurfaceArrayCreator` | `LayerBlueprintNode` binning | Surface lookup acceleration |
| Manual volume construction | `StaticBlueprintNode` | Fixed geometry volumes |
| Manual portal management | Automatic portal system | Inter-volume navigation |

### Include File Updates

```cpp
// Gen1 includes
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"

// Gen3 includes
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
```

## Step-by-Step Migration Examples

### Example 1: Simple Layer Creation

**Gen1 Code:**
```cpp
void createPixelLayerGen1(const GeometryContext& gctx) {
    // Configure layer creator
    Acts::LayerCreator::Config layerConfig;
    layerConfig.surfaceArrayCreator = std::make_shared<Acts::SurfaceArrayCreator>();
    layerConfig.cylinderZtolerance = 5_mm;
    layerConfig.cylinderPhiTolerance = 0.1;
    Acts::LayerCreator layerCreator(layerConfig, logger);
    
    // Create surfaces (implementation not shown)
    auto surfaces = createPixelSurfaces();
    
    // Create layer
    auto layer = layerCreator.cylinderLayer(gctx, surfaces, 16, 20);
    
    // Manual integration into volume (complex process)
    integateLayerIntoVolume(layer);
}
```

**Gen3 Code:**
```cpp
void createPixelLayerGen3(CylinderContainerBlueprintNode& parent) {
    // Create surfaces (same as Gen1)
    auto surfaces = createPixelSurfaces();
    
    // Add layer to parent container
    auto& layer = parent.addLayer("PixelLayer", [&](LayerBlueprintNode& layerNode) {
        layerNode.setLayerType(LayerBlueprintNode::LayerType::Cylinder)
                 .setSurfaces(surfaces)
                 .setBinning(Acts::BinningType::equidistant, Acts::BinningValue::binPhi, 16)
                 .setBinning(Acts::BinningType::equidistant, Acts::BinningValue::binZ, 20)
                 .setEnvelope(Acts::ExtentEnvelope{{.z = {2_mm, 2_mm}}});
    });
    // Automatic integration - no manual work needed!
}
```

### Example 2: Volume Hierarchy Construction

**Gen1 Code:**
```cpp
std::shared_ptr<Acts::TrackingVolume> buildDetectorGen1(const GeometryContext& gctx) {
    // Create helpers with complex configuration
    Acts::CylinderVolumeHelper::Config volConfig;
    volConfig.layerArrayCreator = std::make_shared<Acts::LayerArrayCreator>();
    volConfig.trackingVolumeArrayCreator = std::make_shared<Acts::TrackingVolumeArrayCreator>();
    volConfig.passiveLayerThickness = 5_mm;
    Acts::CylinderVolumeHelper volumeHelper(volConfig, logger);
    
    // Manually create each layer
    auto pixelLayers = createPixelLayers(gctx);
    auto stripLayers = createStripLayers(gctx);
    
    // Manually calculate bounds for each volume
    auto pixelBounds = calculatePixelBounds(pixelLayers);
    auto stripBounds = calculateStripBounds(stripLayers);
    
    // Create individual volumes
    auto pixelVolume = volumeHelper.createTrackingVolume(
        gctx, pixelLayers, pixelBounds, Transform3::Identity(), "Pixel");
    auto stripVolume = volumeHelper.createTrackingVolume(
        gctx, stripLayers, stripBounds, Transform3::Identity(), "Strip");
    
    // Manually create container
    std::vector<Acts::TrackingVolumePtr> volumes = {pixelVolume, stripVolume};
    auto containerBounds = calculateContainerBounds(volumes);
    
    return volumeHelper.createContainerTrackingVolume(
        gctx, volumes, containerBounds, Transform3::Identity(), "Detector");
}
```

**Gen3 Code:**
```cpp
std::unique_ptr<Acts::TrackingGeometry> buildDetectorGen3(const GeometryContext& gctx) {
    using namespace Acts::Experimental;
    
    // Create root blueprint
    Blueprint::Config cfg;
    cfg.envelope = Acts::ExtentEnvelope{{.r = {5_mm, 10_mm}, .z = {10_mm, 10_mm}}};
    auto root = std::make_unique<Blueprint>(cfg);
    
    // Add detector container (automatic sizing)
    auto& detector = root->addCylinderContainer("Detector", Acts::AxisDirection::AxisZ);
    
    // Add pixel detector
    auto& pixel = detector.addCylinderContainer("Pixel", Acts::AxisDirection::AxisZ);
    addPixelLayers(pixel);  // Simple function call
    
    // Add strip detector  
    auto& strip = detector.addCylinderContainer("Strip", Acts::AxisDirection::AxisZ);
    addStripLayers(strip);  // Simple function call
    
    // Automatic construction - handles all sizing, portals, navigation
    return root->construct({}, gctx);
}

void addPixelLayers(CylinderContainerBlueprintNode& pixel) {
    auto& barrel = pixel.addCylinderContainer("PixelBarrel", Acts::AxisDirection::AxisR);
    // Add layers - automatic sizing and portal management
    for (int i = 0; i < 4; ++i) {
        double r = 50_mm + i * 30_mm;
        auto& layer = barrel.addLayer("Layer" + std::to_string(i), [&](auto& ln) {
            ln.setLayerType(LayerBlueprintNode::LayerType::Cylinder)
              .setSurfaces(createPixelSurfaces(r));
        });
    }
}
```

### Example 3: Complex Multi-Subsystem Detector

**Gen1 Code (Simplified):**
```cpp
// Gen1 requires extensive manual management
std::shared_ptr<Acts::TrackingVolume> buildComplexDetectorGen1() {
    // 100+ lines of builder configuration
    // Manual layer creation and bounds calculation
    // Complex portal management between subsystems
    // Error-prone volume assembly
    // Manual navigation policy setup
    
    // This is typically 200-500 lines of complex code
    return finalVolume;
}
```

**Gen3 Code:**
```cpp
std::unique_ptr<Acts::TrackingGeometry> buildComplexDetectorGen3() {
    using namespace Acts::Experimental;
    
    Blueprint::Config cfg;
    cfg.envelope = Acts::ExtentEnvelope{{.r = {10_mm, 20_mm}, .z = {50_mm, 50_mm}}};
    auto root = std::make_unique<Blueprint>(cfg);
    
    // Inner detector
    auto& innerDet = root->addCylinderContainer("InnerDetector", Acts::AxisDirection::AxisZ);
    
    // Pixel system with material designation
    auto& pixelMaterial = innerDet.addMaterial("PixelSilicon", [&](auto& material) {
        auto& pixel = material.addCylinderContainer("Pixel", Acts::AxisDirection::AxisZ);
        buildPixelSystem(pixel);
    });
    
    // Strip system
    auto& strip = innerDet.addCylinderContainer("Strip", Acts::AxisDirection::AxisZ);
    buildStripSystem(strip);
    
    // TRT system with geometry identification
    auto& trt = innerDet.withGeometryIdentifier([&](auto& geoId) {
        auto& trtContainer = geoId.addCylinderContainer("TRT", Acts::AxisDirection::AxisZ);
        buildTRTSystem(trtContainer);
    });
    
    // Calorimeter (cuboid system)
    auto& calo = root->addCuboidContainer("Calorimeter", Acts::AxisDirection::AxisX);
    buildCalorimeterSystem(calo);
    
    // Automatic construction - handles everything!
    return root->construct({}, Acts::GeometryContext{});
}
```

## Common Migration Patterns

### Pattern 1: Layer Creator Replacement

**Before (Gen1):**
```cpp
Acts::LayerCreator::Config config;
config.surfaceArrayCreator = surfaceArrayCreator;
config.cylinderZtolerance = 5_mm;
Acts::LayerCreator creator(config, logger);

auto layer = creator.cylinderLayer(gctx, surfaces, phiBins, zBins);
```

**After (Gen3):**
```cpp
parent.addLayer("MyLayer", [&](LayerBlueprintNode& layer) {
    layer.setLayerType(LayerBlueprintNode::LayerType::Cylinder)
         .setSurfaces(surfaces)
         .setBinning(Acts::BinningType::equidistant, Acts::BinningValue::binPhi, phiBins)
         .setBinning(Acts::BinningType::equidistant, Acts::BinningValue::binZ, zBins)
         .setEnvelope(Acts::ExtentEnvelope{{.z = {2.5_mm, 2.5_mm}}});
});
```

### Pattern 2: Volume Helper Replacement

**Before (Gen1):**
```cpp
Acts::CylinderVolumeHelper::Config config;
config.layerArrayCreator = layerArrayCreator;
config.trackingVolumeArrayCreator = volumeArrayCreator;
Acts::CylinderVolumeHelper helper(config, logger);

auto volume = helper.createTrackingVolume(gctx, layers, bounds, transform, name);
```

**After (Gen3):**
```cpp
auto& container = parent.addCylinderContainer(name, Acts::AxisDirection::AxisZ);
// Add layers to container - automatic volume creation and sizing
for (const auto& layerInfo : layerInfos) {
    container.addLayer(layerInfo.name, [&](auto& layer) {
        layer.setSurfaces(layerInfo.surfaces)
             .setLayerType(layerInfo.type);
    });
}
```

### Pattern 3: Container Volume Creation

**Before (Gen1):**
```cpp
std::vector<Acts::TrackingVolumePtr> childVolumes;
// Manual assembly of child volumes...

auto containerBounds = calculateContainerBounds(childVolumes);
auto container = helper.createContainerTrackingVolume(
    gctx, childVolumes, containerBounds, transform, name);
```

**After (Gen3):**
```cpp
auto& container = parent.addCylinderContainer(name, direction);
container.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::Gap);

// Add children - automatic container sizing
for (const auto& childInfo : childInfos) {
    addChildToContainer(container, childInfo);
}
// No manual bounds calculation or assembly needed
```

## Configuration Migration

### Volume Attachment Strategies

Map your Gen1 volume arrangement logic to Gen3 attachment strategies:

```cpp
// Gen1: Manual volume positioning
auto vol1Transform = Transform3::Identity() * Translation3(Vector3(0, 0, -100_mm));
auto vol2Transform = Transform3::Identity() * Translation3(Vector3(0, 0, 100_mm));

// Gen3: Automatic positioning with strategies
container.setAttachmentStrategy(Acts::CylinderVolumeStack::AttachmentStrategy::Gap);
container.setResizeStrategy(Acts::CylinderVolumeStack::ResizeStrategy::Expand);
```

### Surface Array Configuration

Map your Gen1 surface array configuration to Gen3 layer binning:

```cpp
// Gen1: Surface array creator configuration
Acts::SurfaceArrayCreator::Config arrayConfig;
// Complex configuration...
auto arrayCreator = std::make_shared<Acts::SurfaceArrayCreator>(arrayConfig);

// Gen3: Simple binning configuration
layer.setBinning(Acts::BinningType::equidistant, Acts::BinningValue::binPhi, 16)
     .setBinning(Acts::BinningType::arbitrary, Acts::BinningValue::binR, radialBoundaries);
```

## Testing Your Migration

### 1. Geometry Validation

```cpp
// Validate that Gen3 produces equivalent geometry
auto gen1Geometry = buildDetectorGen1(gctx);
auto gen3Geometry = buildDetectorGen3(gctx);

// Check volume hierarchy
validateVolumeHierarchy(gen1Geometry, gen3Geometry);

// Check portal connectivity  
validatePortalConnectivity(gen3Geometry);

// Check navigation consistency
validateNavigation(gen1Geometry, gen3Geometry);
```

### 2. Performance Comparison

```cpp
// Compare construction performance
auto start = std::chrono::high_resolution_clock::now();
auto gen1Geo = buildDetectorGen1(gctx);
auto gen1Time = std::chrono::high_resolution_clock::now() - start;

start = std::chrono::high_resolution_clock::now();
auto gen3Geo = buildDetectorGen3(gctx);
auto gen3Time = std::chrono::high_resolution_clock::now() - start;

std::cout << "Gen1 construction time: " << gen1Time.count() << "ms\n";
std::cout << "Gen3 construction time: " << gen3Time.count() << "ms\n";
```

### 3. Debugging Tools

```cpp
// Export tree structure for debugging
std::ofstream dotFile("geometry_tree.dot");
blueprint->graphviz(dotFile);

// Use graphviz to visualize: dot -Tpng geometry_tree.dot -o tree.png

// Export geometry for visualization
Acts::ObjVisualization3D vis;
trackingGeometry->visualize(vis, gctx);
std::ofstream objFile("geometry.obj");
vis.write(objFile);
```

## Common Migration Issues

### Issue 1: Missing Portal Connections

**Problem:** Volumes not properly connected in navigation.

**Solution:** Ensure proper container hierarchy and check attachment strategies.

```cpp
// Make sure containers have the right axis direction
auto& container = root->addCylinderContainer("Detector", Acts::AxisDirection::AxisZ);
//                                                        ^^^^^^^^^^^^^^^^^^^^^^^^
//                                                        Critical for proper stacking
```

### Issue 2: Incorrect Volume Sizing

**Problem:** Volumes are too small or too large.

**Solution:** Check envelope configuration and resize strategies.

```cpp
// Configure appropriate envelopes
Blueprint::Config cfg;
cfg.envelope = Acts::ExtentEnvelope{{
    .r = {clearance_inner, clearance_outer},
    .z = {clearance_neg, clearance_pos}
}};
```

### Issue 3: Layer Surface Integration

**Problem:** Surfaces not properly integrated into layers.

**Solution:** Ensure proper layer type and surface assignment.

```cpp
layer.setLayerType(LayerBlueprintNode::LayerType::Cylinder)  // Match geometry type
     .setSurfaces(surfaces)                                 // Must be non-empty
     .setTransform(appropriateTransform);                   // Must match surface positions
```

## Benefits After Migration

After successful migration to Gen3, you'll experience:

1. **Reduced Code Complexity**: 50-80% reduction in geometry construction code
2. **Better Maintainability**: Clear hierarchical structure
3. **Automatic Portal Management**: No more manual portal setup
4. **Intelligent Sizing**: Automatic bounds calculation
5. **Better Debugging**: Tree visualization and clear error messages
6. **Improved Performance**: Optimized construction process
7. **Future-Proof**: Foundation for upcoming geometry features

The migration effort is typically recouped within weeks through reduced debugging time and easier geometry modifications.