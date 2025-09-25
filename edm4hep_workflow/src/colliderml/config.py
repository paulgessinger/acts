from enum import StrEnum
from pathlib import Path
from typing import Any, Self
import pydantic
import tomllib


class TomlConfigBase(pydantic.BaseModel):
    @classmethod
    def load(cls, filename: Path | None) -> Self:
        if filename is None:
            return cls()

        assert filename.exists(), f"Config file {filename} does not exist"

        with filename.open("rb") as f:
            raw = tomllib.load(f)

        return cls.model_validate(raw)


class SimHitReading(pydantic.BaseModel):
    particleR: tuple[float | None, float | None] = (None, None)
    particleZ: tuple[float | None, float | None] = (None, None)


class Config(TomlConfigBase):

    sim_hit_reading: SimHitReading = pydantic.Field(default_factory=SimHitReading)


class CardCustomizations(pydantic.BaseModel):
    run_card: dict[str, Any] = pydantic.Field(default_factory=dict)
    shower_card: dict[str, Any] = pydantic.Field(default_factory=dict)
    pythia8_card: dict[str, Any] = pydantic.Field(default_factory=dict)


class RunMode(StrEnum):
    nlo_fxfx = "nlo_fxfx"
    lo_mlm = "lo_mlm"


class SampleConfig(TomlConfigBase):
    label: str

    model: str
    run_mode: RunMode

    definitions: str = ""

    generate_command: str

    card_customizations: CardCustomizations = pydantic.Field(
        default_factory=CardCustomizations
    )
