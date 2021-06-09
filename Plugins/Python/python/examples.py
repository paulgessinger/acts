from acts._acts._examples import *
import acts._acts._examples

from acts import _acts

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