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


class BaseModel(pydantic.BaseModel):
    class Config:
        extra = "forbid"


class BinningValue(Enum):
    binX = "binX"
    binY = "binY"
    binZ = "binZ"
    binR = "binR"
    binPhi = "binPhi"
    binRPhi = "binRPhi"
    binH = "binH"
    binEta = "binEta"
    binMag = "binMag"


class BinningOption(Enum):
    open = "open"
    closed = "closed"


class BinningType(Enum):
    equidistant = "equidistant"
    arbitrary = "arbitrary"


class BinDataSingle(BaseModel):
    min: float
    max: float
    bins: int
    value: BinningValue
    option: BinningOption
    type: BinningType


class BinDataWithSub(BinDataSingle):
    subdata: Optional[Union["BinDataSingle", "BinDataWithSub"]]
    subadditive: Optional[bool]


BinDataSingle.update_forward_refs()
BinDataWithSub.update_forward_refs()

BinData = Union[BinDataSingle, BinDataWithSub]
