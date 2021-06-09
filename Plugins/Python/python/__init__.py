from acts._version import __version__

from acts._acts import *
import acts._acts as _acts

def Propagator(stepper, nav):
    for prefix in ("Eigen", "Atlas", "StraightLine"):
        _stepper = getattr(_acts, f"{prefix}Stepper")
        if isinstance(stepper, _stepper):
            return getattr(_acts._propagator, f"{prefix}Propagator")(stepper, nav)
    raise TypeError(f"Unknown stepper {type(stepper).__name__}")
