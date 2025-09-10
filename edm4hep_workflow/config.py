from pathlib import Path
import pydantic
import tomllib


class SimHitReading(pydantic.BaseModel):
    particleR: tuple[float | None, float | None] = (None, None)
    particleZ: tuple[float | None, float | None] = (None, None)


class Config(pydantic.BaseModel):
    sim_hit_reading: SimHitReading = pydantic.Field(default_factory=SimHitReading)

    @staticmethod
    def load(filename: Path | None) -> "Config":
        if filename is None:
            return Config()

        assert filename.exists(), f"Config file {filename} does not exist"

        with filename.open("rb") as f:
            raw = tomllib.load(f)

        return Config.model_validate(raw)
