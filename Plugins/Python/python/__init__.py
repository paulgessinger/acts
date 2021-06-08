from acts._version import __version__

from acts._acts import *
import acts._acts as _acts

import pkgutil
import importlib
# for importer, modname, ispkg in pkgutil.iter_modules(__path__, __name__+"."):
    # print(importer, modname, ispkg)
    # importlib.import_module()

def Propagator(stepper, nav):
    if isinstance(stepper, _acts.EigenStepper):
        return _acts.propagator.EigenPropagator(stepper, nav)
    elif isinstance(stepper, _acts.AtlasStepper):
        return _acts.propagator.AtlasPropagator(stepper, nav)
    elif isinstance(stepper, _acts.StraightLineStepper):
        return _acts.propagator.StraightLinePropagator(stepper, nav)
    else:
        raise TypeError(f"Unknown stepper {type(stepper).__name__}")

class _propagationAlgorithmHelper():

    @staticmethod
    def Config(propagator):
        if isinstance(propagator, _acts.propagator.EigenPropagator):
            return _acts.propagator.EigenPropagationAlgorithm.Config(propagator)
        elif isinstance(propagator, _acts.propagator.AtlasPropagator):
            return _acts.propagator.AtlasPropagationAlgorithm.Config(propagator)
        elif isinstance(propagator, _acts.propagator.StraightLinePropagator):
            return _acts.propagator.StraightLinePropagationAlgorithm.Config(propagator)
        else:
            raise TypeError(f"Unknown propagator {type(propagator).__name__}")


    @staticmethod
    def __call__(config, *args):
        if isinstance(config, _acts.propagator.EigenPropagationAlgorithm.Config):
            return _acts.propagator.EigenPropagationAlgorithm(config, *args)
        elif isinstance(config, _acts.propagator.AtlasPropagationAlgorithm.Config):
            return _acts.propagator.AtlasPropagationAlgorithm(config, *args)
        elif isinstance(config, _acts.propagator.StraightLinePropagationAlgorithm.Config):
            return _acts.propagator.StraightLinePropagationAlgorithm(config, *args)
        else:
            raise TypeError(f"Unknown propagator config {type(config).__name__}")

PropagationAlgorithm = _propagationAlgorithmHelper()