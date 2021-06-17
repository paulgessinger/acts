from ._acts._examples import *
from . import _acts
from ._acts import _examples
import acts

import inspect

_propagators = []
_propagation_algorithms = []
for prefix in ("Eigen", "Atlas", "StraightLine"):
    _propagators.append(getattr(_acts._propagator, f"{prefix}Propagator"))
    _propagation_algorithms.append(
        getattr(_acts._examples, f"{prefix}PropagationAlgorithm")
    )


class _propagationAlgorithmHelper:
    @staticmethod
    def Config(propagator, **kwargs):
        for prop, alg in zip(_propagators, _propagation_algorithms):
            if isinstance(propagator, prop):
                cfg = alg.Config(propagator)
                for k, v in kwargs.items():
                    setattr(cfg, k, v)
                return cfg

        raise TypeError(f"Unknown propagator {type(propagator).__name__}")

    @staticmethod
    def __call__(*args, **kwargs):
        config_or_propagator = kwargs.get("propagator")
        if config_or_propagator is None:
            config_or_propagator = args[0]
            args = args[1:]
        for prop, alg in zip(_propagators, _propagation_algorithms):
            if isinstance(config_or_propagator, alg.Config):
                return alg(config_or_propagator, *args, **kwargs)
            if isinstance(config_or_propagator, prop):
                cfg = alg.Config(config_or_propagator)
                _kwargs = {}
                for k, v in kwargs.items():
                    if hasattr(cfg, k):
                        setattr(cfg, k, v)
                    else:
                        _kwargs[k] = v
                return alg(cfg, *args, **_kwargs)

        raise TypeError(f"Unknown propagator config {type(config).__name__}")


PropagationAlgorithm = _propagationAlgorithmHelper()


def _detector_create(cls):
    def create(*args, mdecorator=None, **kwargs):
        cfg = cls.Config()
        _kwargs = {}
        for k, v in kwargs.items():
            if hasattr(cfg, k):
                setattr(cfg, k, v)
            else:
                _kwargs[k] = v
        det = cls()
        tg, deco = det.finalize(cfg, mdecorator, *args, **_kwargs)
        return det, tg, deco

    return create


for m in (_acts._examples, _acts._examples.dd4hep):
    for name, cls in inspect.getmembers(m, inspect.isclass):
        # print(name)
        if not hasattr(cls, "Config"):
            continue
        if name.endswith("Detector"):
            cls.create = _detector_create(cls)
        else:
            cls.__init__ = acts._make_config_adapter(cls.__init__)
