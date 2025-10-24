from enum import StrEnum
from pathlib import Path
from typing import Annotated, Any, Literal, Self
import pydantic
from pydantic import BeforeValidator
import tomllib
import acts


def _parse_pdg_particle(value: str | int) -> int:
    """
    Parse a PDG particle from either a string name or an integer PDG code.

    Examples:
        "mu-" -> 13
        "e+" -> -11
        2212 -> 2212
    """
    if isinstance(value, int):
        return value

    if isinstance(value, str):

        pdg_particle = acts.parsePdgParticle(value)
        return int(pdg_particle)

    raise ValueError(
        f"PDG particle must be a string name or integer, got {type(value)}"
    )


# Type annotation for PDG particle fields that accepts both string names and integers
PdgParticle = Annotated[int, BeforeValidator(_parse_pdg_particle)]


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


class ParticleHandlerType(StrEnum):
    default = "default"
    full_truth = "Geant4FullTruthParticleHandler"
    tracking_cylinder = "Geant4TCUserParticleHandler"


class SimulationConfig(pydantic.BaseModel):
    minimal_kinetic_energy: float = 0.1  # GeV

    user_particle_handler: ParticleHandlerType = ParticleHandlerType.full_truth


class DigitizationConfig(pydantic.BaseModel):
    config_file: str


class Config(TomlConfigBase):
    sim_hit_reading: SimHitReading = pydantic.Field(default_factory=SimHitReading)

    simulation: SimulationConfig = pydantic.Field(default_factory=SimulationConfig)

    digitization: DigitizationConfig = pydantic.Field(
        default_factory=DigitizationConfig
    )


class CardCustomizations(pydantic.BaseModel):
    run_card: dict[str, str] = pydantic.Field(default_factory=dict)
    shower_card: dict[str, str] = pydantic.Field(default_factory=dict)
    pythia8_card: dict[str, str] = pydantic.Field(default_factory=dict)


class RunMode(StrEnum):
    nlo_fxfx = "nlo_fxfx"
    lo_mlm = "lo_mlm"


class SampleConfigBase(TomlConfigBase):
    label: str

    ebeam1: float
    ebeam2: float


class MadgraphSampleConfig(SampleConfigBase):

    model: str
    run_mode: RunMode

    definitions: str = ""

    generate_command: str

    # Speculatively produce more events to counter vetoes etc.
    nevents_scale_factor: float = 1.0

    card_customizations: CardCustomizations = pydantic.Field(
        default_factory=CardCustomizations
    )

    @property
    def run_card(self) -> dict[str, str]:
        return {
            **self.card_customizations.run_card,
            "ebeam1": str(self.ebeam1),
            "ebeam2": str(self.ebeam2),
        }


class Pythia8SampleConfig(SampleConfigBase):
    settings: list[str]

    pdg_beam0_internal: PdgParticle = pydantic.Field(alias="pdg_beam0")
    pdg_beam1_internal: PdgParticle = pydantic.Field(alias="pdg_beam1")

    @property
    def pdg_beam0(self) -> acts.PdgParticle:
        return acts.PdgParticle(self.pdg_beam0_internal)

    @property
    def pdg_beam1(self) -> acts.PdgParticle:
        return acts.PdgParticle(self.pdg_beam1_internal)

    @property
    def cms_energy(self) -> float:
        return self.ebeam1 + self.ebeam2


class PileupStrategy(StrEnum):
    fixed = "fixed"
    poisson = "poisson"


Unit = Literal["mm", "um", "m", "cm", "ns", "s", "ps"]


class ValueWithUnit(pydantic.BaseModel):
    value: float
    unit: Unit

    def to_acts(self) -> float:
        return self.value * getattr(acts.UnitConstants, self.unit)

    def __str__(self) -> str:
        return f"{self.value} {self.unit}"


class Vector4WithUnit(
    pydantic.RootModel[
        tuple[ValueWithUnit, ValueWithUnit, ValueWithUnit, ValueWithUnit]
    ]
):

    def to_acts(self) -> acts.Vector4:
        return acts.Vector4(*map(lambda v: v.to_acts(), self.root))

    @staticmethod
    def zeros():
        return Vector4WithUnit(
            ValueWithUnit(value=0.0, unit="mm"),
            ValueWithUnit(value=0.0, unit="mm"),
            ValueWithUnit(value=0.0, unit="mm"),
            ValueWithUnit(value=0.0, unit="ns"),
        )

    def __str__(self) -> str:
        x, y, z, t = self.root
        return f"[{x}, {y}, {z}, {t}]"


class PileupConfig(TomlConfigBase):
    # vertex_units: tuple[str, str, str, str] = ("mm", "mm", "mm", "ns")
    # vertex_mean: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
    vertex_mean: Vector4WithUnit = pydantic.Field(default_factory=Vector4WithUnit.zeros)
    vertex_stddev: Vector4WithUnit = pydantic.Field(
        default_factory=Vector4WithUnit.zeros
    )

    strategy: PileupStrategy = PileupStrategy.poisson
