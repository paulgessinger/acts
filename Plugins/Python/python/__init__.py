from acts._version import __version__

from acts._acts import *
import acts._acts as _acts

import acts.examples

def Propagator(stepper, nav):
    if isinstance(stepper, _acts.EigenStepper):
        return _acts.propagator.EigenPropagator(stepper, nav)
    elif isinstance(stepper, _acts.AtlasStepper):
        return _acts.propagator.AtlasPropagator(stepper, nav)
    elif isinstance(stepper, _acts.StraightLineStepper):
        return _acts.propagator.StraightLinePropagator(stepper, nav)
    else:
        raise TypeError(f"Unknown stepper {type(stepper).__name__}")
