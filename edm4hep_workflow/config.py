from pathlib import Path
from typing import Self
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


class CardCustomization(pydantic.BaseModel):
    pass


class CardCustomizations(pydantic.BaseModel):
    run_card: CardCustomization = pydantic.Field(default_factory=CardCustomization)
    shower_card: CardCustomization = pydantic.Field(default_factory=CardCustomization)
    pythia8_card: CardCustomization = pydantic.Field(default_factory=CardCustomization)


class SampleConfig(TomlConfigBase):

    model: str

    definitions: list[str] | None = None

    generate_command: str
