from .ActsPythonBindings import *
from .ActsPythonBindings import __version__
from . import ActsPythonBindings
from ._adapter import _patch_config

print(__file__)


def Propagator(stepper, navigator):
    for prefix in ("Eigen", "Atlas", "StraightLine"):
        _stepper = getattr(ActsPythonBindings, f"{prefix}Stepper")
        if isinstance(stepper, _stepper):
            return getattr(ActsPythonBindings._propagator, f"{prefix}Propagator")(
                stepper, navigator
            )
    raise TypeError(f"Unknown stepper {type(stepper).__name__}")


_patch_config(ActsPythonBindings)
