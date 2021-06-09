from acts._acts._examples import *
import acts._acts._examples

from acts import _acts

_propagators = []
_propagation_algorithms = []
for prefix in ("Eigen", "Atlas", "StraightLine"):
    _propagators.append(getattr(acts._acts._propagator, f"{prefix}Propagator"))
    _propagation_algorithms.append(getattr(acts._acts._examples, f"{prefix}PropagationAlgorithm"))


class _propagationAlgorithmHelper():

    @staticmethod
    def Config(propagator):
        for prop, alg in zip(_propagators, _propagation_algorithms):
            if isinstance(propagator, prop):
                return alg.Config(propagator)

        raise TypeError(f"Unknown propagator {type(propagator).__name__}")


    @staticmethod
    def __call__(config, *args):
        for prop, alg in zip(_propagators, _propagation_algorithms):
            if isinstance(config, alg.Config):
                return alg(config, *args)

        raise TypeError(f"Unknown propagator config {type(config).__name__}")

PropagationAlgorithm = _propagationAlgorithmHelper()