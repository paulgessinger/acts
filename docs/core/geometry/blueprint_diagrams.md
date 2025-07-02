# Blueprint System Diagrams

This document provides visual explanations of the Blueprint system's architecture, construction phases, and data flow.

## Blueprint Tree Structure

### Hierarchical Organization

```
                           Blueprint (Root)
                          /       |        \
                    Container   Material   GeometryID
                    Nodes       Designator  Node
                   /     \      Node        |
            Cylinder    Cuboid     |        └── Child Nodes
            Container   Container  |            (inherit IDs)
            /     \     /     \    |
        Static  Layer Static Layer └── Child Nodes
        Volume  Node  Volume Node      (inherit material)
```

### Example: Complex Detector Hierarchy

```
Acts Tracking Geometry
└── Blueprint (Root)
    └── InnerDetector (CylinderContainer, AxisZ)
        ├── Pixel (CylinderContainer, AxisZ)
        │   ├── PixelBarrel (CylinderContainer, AxisR)
        │   │   ├── PixelBarrelLayer0 (LayerNode)
        │   │   ├── PixelBarrelLayer1 (LayerNode)
        │   │   ├── PixelBarrelLayer2 (LayerNode)
        │   │   └── PixelBarrelLayer3 (LayerNode)
        │   ├── PixelPosEndcap (CylinderContainer, AxisZ)
        │   │   ├── PixelPosEndcapDisk0 (LayerNode)
        │   │   ├── PixelPosEndcapDisk1 (LayerNode)
        │   │   └── PixelPosEndcapDisk2 (LayerNode)
        │   └── PixelNegEndcap (CylinderContainer, AxisZ)
        │       ├── PixelNegEndcapDisk0 (LayerNode)
        │       ├── PixelNegEndcapDisk1 (LayerNode)
        │       └── PixelNegEndcapDisk2 (LayerNode)
        ├── Strip (CylinderContainer, AxisZ)
        │   ├── StripBarrel (CylinderContainer, AxisR)
        │   │   ├── StripBarrelLayer0 (LayerNode)
        │   │   ├── StripBarrelLayer1 (LayerNode)
        │   │   ├── StripBarrelLayer2 (LayerNode)
        │   │   └── StripBarrelLayer3 (LayerNode)
        │   ├── StripPosEndcap (CylinderContainer, AxisZ)
        │   └── StripNegEndcap (CylinderContainer, AxisZ)
        └── TRT_MaterialRegion (MaterialDesignatorNode)
            └── TRT (CylinderContainer, AxisZ)
                ├── TRTBarrel (StaticVolume)
                ├── TRTPosEndcap (StaticVolume)
                └── TRTNegEndcap (StaticVolume)
```

## Three-Phase Construction Process

### Phase 1: Build - Volume Sizing and Bounds Calculation

```
Phase 1: BUILD
==============

Step 1: Leaf nodes compute their volumes
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│ LayerNode       │    │ StaticVolume    │    │ LayerNode       │
│ "PixelLayer0"   │    │ "TRTBarrel"     │    │ "StripLayer0"   │
│ ┌─────────────┐ │    │ ┌─────────────┐ │    │ ┌─────────────┐ │
│ │Volume:      │ │    │ │Volume:      │ │    │ │Volume:      │ │
│ │r: 50-55mm   │ │    │ │r: 600-900mm │ │    │ │r: 200-205mm │ │
│ │z: ±300mm    │ │    │ │z: ±800mm    │ │    │ │z: ±400mm    │ │
│ └─────────────┘ │    │ └─────────────┘ │    │ └─────────────┘ │
└─────────────────┘    └─────────────────┘    └─────────────────┘

Step 2: Container nodes compute bounding volumes
┌─────────────────────────────────────────────────────────────────┐
│ CylinderContainer "PixelBarrel" (AxisR)                        │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ Computed Bounds:                                            │ │
│ │ r: 48-157mm (envelope + child extents)                     │ │
│ │ z: ±302mm (envelope + child extents)                       │ │
│ │                                                             │ │
│ │ Child Volumes:                                              │ │
│ │ ├── PixelLayer0: r: 50-55mm, z: ±300mm                     │ │
│ │ ├── PixelLayer1: r: 80-85mm, z: ±300mm                     │ │
│ │ ├── PixelLayer2: r: 110-115mm, z: ±300mm                   │ │
│ │ └── PixelLayer3: r: 140-145mm, z: ±300mm                   │ │
│ └─────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘

Step 3: Size propagation up the tree
┌─────────────────────────────────────────────────────────────────┐
│ Blueprint (Root)                                                │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ World Volume Bounds:                                        │ │
│ │ r: 0-920mm (largest child + envelope)                      │ │
│ │ z: ±820mm (largest child + envelope)                       │ │
│ └─────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
```

### Phase 2: Connect - Portal Creation and Inter-Volume Connections

```
Phase 2: CONNECT
================

Step 1: Create boundary surfaces for each volume
┌─────────────────┐         ┌─────────────────┐
│ Volume A        │         │ Volume B        │
│ ┌─────────────┐ │         │ ┌─────────────┐ │
│ │    ┌───┐    │ │         │ │    ┌───┐    │ │
│ │    │   │ ◄──┼─┼─────────┼─┼──► │   │    │ │
│ │    └───┘    │ │ Shared  │ │    └───┘    │ │
│ └─────────────┘ │ Boundary│ └─────────────┘ │
└─────────────────┘         └─────────────────┘

Step 2: Convert boundaries to portals with links
                    Portal
                   ┌─────┐
     Volume A ◄────┤ │ │ ├────► Volume B
                   └─────┘
                  TrivialPortalLink
                  (one-to-one)

Step 3: Create portal shells for complex arrangements
Container Node Portal Shell:
┌─────────────────────────────────────────────────────────────────┐
│ CylinderPortalShell                                             │
│                                                                 │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │                    Outer Cylinder                          │ │
│ │  ┌─────────────────────────────────────────────────────┐   │ │
│ │  │                Inner Cylinder                       │   │ │
│ │  │   ┌─────┐  ┌─────┐  ┌─────┐  ┌─────┐              │   │ │
│ │  │   │ V1  │  │ V2  │  │ V3  │  │ V4  │              │   │ │
│ │  │   └─────┘  └─────┘  └─────┘  └─────┘              │   │ │
│ │  └─────────────────────────────────────────────────────┘   │ │
│ │              Negative Z    Positive Z                      │ │
│ │               Disc         Disc                            │ │
│ └─────────────────────────────────────────────────────────────┘ │
│                                                                 │
│ Portals: Inner Cylinder, Outer Cylinder, Negative Z, Positive Z│
└─────────────────────────────────────────────────────────────────┘
```

### Phase 3: Finalize - TrackingGeometry Assembly

```
Phase 3: FINALIZE
=================

Step 1: Register portals with their volumes
┌─────────────────┐
│ TrackingVolume  │
│ ┌─────────────┐ │    ┌─────────────────┐
│ │             │ │    │ Navigation      │
│ │   Volume    │ │◄───┤ Policy          │
│ │   Geometry  │ │    │ (Surface Array, │
│ │             │ │    │  Try All, etc.) │
│ └─────────────┘ │    └─────────────────┘
│ ┌─────────────┐ │
│ │ Portal      │ │
│ │ Registry    │ │
│ └─────────────┘ │
└─────────────────┘

Step 2: Create navigation policies
Surface Array Policy:
┌─────────────────────────────────────────────┐
│ Surface Grid (for efficient lookup)         │
│ ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┐ │
│ │ S1  │ S2  │ S3  │ S4  │ S5  │ S6  │ S7  │ │
│ ├─────┼─────┼─────┼─────┼─────┼─────┼─────┤ │
│ │ S8  │ S9  │ S10 │ S11 │ S12 │ S13 │ S14 │ │
│ ├─────┼─────┼─────┼─────┼─────┼─────┼─────┤ │
│ │ S15 │ S16 │ S17 │ S18 │ S19 │ S20 │ S21 │ │
│ └─────┴─────┴─────┴─────┴─────┴─────┴─────┘ │
└─────────────────────────────────────────────┘

Step 3: Assemble final TrackingGeometry
┌─────────────────────────────────────────────────────────────────┐
│ TrackingGeometry (Gen3)                                         │
│ ┌─────────────────────────────────────────────────────────────┐ │
│ │ Highest Tracking Volume (World)                             │ │
│ │ ┌─────────────┬─────────────┬─────────────┬─────────────┐   │ │
│ │ │ Volume 1    │ Volume 2    │ Volume 3    │ Volume 4    │   │ │
│ │ │ (Pixel)     │ (Strip)     │ (TRT)       │ (Calo)      │   │ │
│ │ │ ┌─────────┐ │ ┌─────────┐ │ ┌─────────┐ │ ┌─────────┐ │   │ │
│ │ │ │Surfaces │ │ │Surfaces │ │ │Material │ │ │Material │ │   │ │
│ │ │ │Navigation│ │ │Navigation│ │ │Only     │ │ │Only     │ │   │ │
│ │ │ └─────────┘ │ └─────────┘ │ └─────────┘ │ └─────────┘ │   │ │
│ │ └─────────────┴─────────────┴─────────────┴─────────────┘   │ │
│ └─────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
```

## Volume Stacking Strategies

### Cylinder Container Stacking (AxisZ)

```
Attachment Strategy: Gap
========================

Input Volumes:
┌─────┐     ┌─────┐     ┌─────┐
│ V1  │     │ V2  │     │ V3  │
│z:-50│     │z:100│     │z:300│
│  to │     │ to │     │ to │
│ +50 │     │200 │     │400 │
└─────┘     └─────┘     └─────┘

Result with Gap Strategy:
┌─────┐┌───┐┌─────┐┌───┐┌─────┐
│ V1  ││Gap││ V2  ││Gap││ V3  │
│     ││ 1 ││     ││ 2 ││     │
└─────┘└───┘└─────┘└───┘└─────┘
-100  -50  0  100 200 250 400 450
      ^     ^     ^     ^
   Gaps inserted automatically

Attachment Strategy: First
==========================

Result with First Strategy:
┌─────┐┌─────┐┌─────┐
│ V1  ││ V2  ││ V3  │
│     ││     ││     │
└─────┘└─────┘└─────┘
-50  +50   150   250
     ^      ^     ^
   All volumes attached to first
```

### Cylinder Container Stacking (AxisR)

```
Radial Stacking (Barrel Layers):
=================================

Input:
Layer 1: r=50-55mm
Layer 2: r=80-85mm  
Layer 3: r=110-115mm
Layer 4: r=140-145mm

Result:
              ┌─────────────────────────┐ ← Layer 4 (r=140-145)
              │┌───────────────────────┐│
              ││┌─────────────────────┐││ ← Layer 3 (r=110-115)
              │││┌───────────────────┐│││
              ││││┌─────────────────┐││││ ← Layer 2 (r=80-85)
              │││││┌───────────────┐│││││
              ││││││               ││││││ ← Layer 1 (r=50-55)
              ││││││     Beam      ││││││
              ││││││     Pipe      ││││││
              │││││└───────────────┘│││││
              ││││└─────────────────┘││││
              │││└───────────────────┘│││
              ││└─────────────────────┘││
              │└───────────────────────┘│
              └─────────────────────────┘

Automatic radial envelope calculation ensures proper spacing.
```

## Portal Connection Patterns

### Simple Portal Connection

```
Volume A                    Volume B
┌──────────────┐           ┌──────────────┐
│              │           │              │
│              ├─────────→ │              │
│              │  Portal   │              │
│              │           │              │
└──────────────┘           └──────────────┘

TrivialPortalLink: Direct one-to-one connection
```

### Grid Portal Connection

```
Container Volume
┌─────────────────────────────────────┐
│  ┌─────┐  ┌─────┐  ┌─────┐  ┌─────┐ │
│  │ V1  │  │ V2  │  │ V3  │  │ V4  │ │
│  └─────┘  └─────┘  └─────┘  └─────┘ │
│  ┌─────┐  ┌─────┐  ┌─────┐  ┌─────┐ │
│  │ V5  │  │ V6  │  │ V7  │  │ V8  │ │
│  └─────┘  └─────┘  └─────┘  └─────┘ │
└─────────────────────────────────────┘
              │
              ▼
GridPortalLink: Uses 2D grid to route between appropriate volumes
┌─────┬─────┬─────┬─────┐
│ V1  │ V2  │ V3  │ V4  │
├─────┼─────┼─────┼─────┤
│ V5  │ V6  │ V7  │ V8  │
└─────┴─────┴─────┴─────┘
```

### Composite Portal Connection

```
Complex Multi-Region Portal:
┌─────────────────────────────────────┐
│             Region A                │
├─────────────────────────────────────┤
│             Region B                │ ── CompositePortalLink
├─────────────────────────────────────┤    (routes based on position)
│             Region C                │
└─────────────────────────────────────┘

Each region connects to different destination volumes based on particle position.
```

## Blueprint System Data Flow

```
User Code Input → Blueprint Tree → Construction Phases → TrackingGeometry

┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│ User constructs │    │ Blueprint       │    │ TrackingGeometry│
│ Blueprint tree  │───→│ Construction    │───→│ (Gen3)          │
│ via API calls   │    │ (3 phases)      │    │ Ready for use   │
└─────────────────┘    └─────────────────┘    └─────────────────┘
         │                       │                       │
         ▼                       ▼                       ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│ • addContainer  │    │ Phase 1: Build  │    │ • Navigation    │
│ • addLayer      │    │ Phase 2: Connect│    │ • Portal system │
│ • addMaterial   │    │ Phase 3:Finalize│    │ • Tracking ready│
│ • addStaticVol  │    │                 │    │                 │
└─────────────────┘    └─────────────────┘    └─────────────────┘
```

These diagrams illustrate how the Blueprint system organizes geometry construction through a clear tree hierarchy, automated three-phase building process, and intelligent portal management system.