#!/usr/bin/env python3

import pydantic
from pydantic import Field

from enum import Enum
import json

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from typing import List, Optional, Tuple, Union

from schema_base import BaseModel, BinData


class Header(BaseModel):
    format_version: Literal[0] = Field(..., alias="format-version")
    value_identifier: Literal["digitization-configuration"] = Field(
        ..., alias="value-identifier"
    )


class Type(Enum):
    Gauss = 1
    GaussTrunk = 2
    GaussClipped = 3
    Uniform = 4
    Digital = 5


class Gauss(BaseModel):
    index: int
    mean: float
    stddev: float
    type: Literal["Gauss"]


class GaussTrunc(BaseModel):
    index: int
    mean: float
    stddev: float
    range: Tuple[float, float]
    type: Literal["GaussTrunc"]


class GaussClipped(BaseModel):
    index: int
    mean: float
    stddev: float
    range: Tuple[float, float]
    max_attempts: int
    type: Literal["GaussClipped"]


class Uniform(BaseModel):
    type: Literal["Uniform"]
    bindata: BinData


class Smearing(BaseModel):
    smearing: List[Union[Gauss, GaussTrunc, GaussClipped, Uniform]]


class Entry(BaseModel):
    class Config:
        extra = "forbid"
        schema_extra = {
            "anyOf": [
                {"required": ["volume"]},
                {"required": ["layer"]},
                {"required": ["sensitive"]},
            ]
        }

    volume: Optional[int] = None
    layer: Optional[int] = None
    sensitive: Optional[int] = None
    value: Smearing


class DigitizationConfig(BaseModel):

    acts_geometry_hierarchy_map: Header = Field(
        ..., alias="acts-geometry-hierarchy-map"
    )
    entries: List[Entry]


# print(DigitizationConfig.parse_file("digi-config-out.json"))

schema = DigitizationConfig.schema()
schema["$schema"] = "https://json-schema.org/draft/2020-12/schema"
print(json.dumps(schema, indent=2))
