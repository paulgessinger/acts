#!/usr/bin/env python3

import pydantic
from pydantic import Field, conlist

from enum import Enum
import json

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from typing import List, Optional, Tuple, Union, Any

from schema_base import BaseModel, BinData


class HeaderVolume(BaseModel):
    format_version: Literal[0] = Field(..., alias="format-version")
    value_identifier: Literal["Material Volume Map"] = Field(
        ..., alias="value-identifier"
    )


class HeaderSurface(BaseModel):
    format_version: Literal[0] = Field(..., alias="format-version")
    value_identifier: Literal["Material Surface Map"] = Field(
        ..., alias="value-identifier"
    )


class Transform(BaseModel):
    rotation: Optional[conlist(float, min_items=9, max_items=9)]
    translation: Optional[conlist(float, min_items=3, max_items=3)]

    class Config:
        extra = "forbid"

        def schema_extra(schema, model):
            schema["properties"]["translation"].update({"type": ["null", "array"]})
            schema["properties"]["rotation"].update({"type": ["null", "array"]})


class BoundsType(Enum):
    ConeBounds = "ConeBounds"
    CylinderBounds = "CylinderBounds"
    DiamondBounds = "DiamondBounds"
    RadialBounds = "RadialBounds"
    EllipseBounds = "EllipseBounds"
    LineBounds = "LineBounds"
    RectangleBounds = "RectangleBounds"
    TrapezoidBounds = "TrapezoidBounds"
    TriangleBounds = "TriangleBounds"
    DiscTrapezoidBounds = "DiscTrapezoidBounds"
    ConvexPolygonBounds = "ConvexPolygonBounds"
    AnnulusBounds = "AnnulusBounds"
    OtherBounds = "OtherBounds"


class Bounds(BaseModel):
    type: BoundsType
    values: List[float]  # this should likely be constrained further


class SurfaceType(Enum):
    ConeSurface = "ConeSurface"
    CylinderSurface = "CylinderSurface"
    DiscSurface = "DiscSurface"
    PerigeeSurface = "PerigeeSurface"
    PlaneSurface = "PlaneSurface"
    StrawSurface = "StrawSurface"
    CurvilinearSurface = "CurvilinearSurface"


class SurfaceValue(BaseModel):
    material: Optional[Any]
    bounds: Bounds
    geo_id: int
    transform: Transform
    type: SurfaceType


class SurfaceEntry(BaseModel):
    class Config:
        extra = "forbid"
        schema_extra = {
            "anyOf": [
                {"required": ["volume"]},
                {"required": ["layer"]},
                {"required": ["sensitive"]},
                {"required": ["boundary"]},
                {"required": ["approach"]},
            ]
        }

    volume: Optional[int] = None
    layer: Optional[int] = None
    sensitive: Optional[int] = None
    boundary: Optional[int] = None
    approach: Optional[int] = None
    value: SurfaceValue


class VolumeBinUtility(BaseModel):
    binningdata: List[BinData]


class VolumeMaterialType(Enum):
    proto = "proto"
    homogeneous = "homogeneous"
    interpolated2D = "interpolated2D"
    interpolated3D = "interpolated3D"


class VolumeMaterial(BaseModel):
    mapMaterial: bool
    binUtility: VolumeBinUtility
    type: VolumeMaterialType


class VolumeValue(BaseModel):
    NAME: str
    material: Optional[VolumeMaterial]


class VolumeEntry(BaseModel):
    value: VolumeValue
    volume: int


class VolumeConfig(BaseModel):
    acts_geometry_hierarchy_map: HeaderVolume = Field(
        ..., alias="acts-geometry-hierarchy-map"
    )
    entries: List[VolumeEntry]


class SurfaceConfig(BaseModel):
    acts_geometry_hierarchy_map: HeaderSurface = Field(
        ..., alias="acts-geometry-hierarchy-map"
    )
    entries: List[SurfaceEntry]


class TopModel(BaseModel):
    Volumes: Any
    Surfaces: SurfaceConfig


# print(
#     TopModel.parse_file("thirdparty/OpenDataDetector/config/odd-material-mapping.json")
# )

schema = TopModel.schema()
schema["$schema"] = "https://json-schema.org/draft/2020-12/schema"
print(json.dumps(schema, indent=2))
