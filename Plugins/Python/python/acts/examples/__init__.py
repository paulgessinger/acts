from acts.ActsPythonBindings._examples import *
from acts import ActsPythonBindings
from acts.ActsPythonBindings import _examples
import acts
from acts._adapter import _patch_config

import inspect

_propagators = []
_propagation_algorithms = []
for prefix in ("Eigen", "Atlas", "StraightLine"):
    _propagators.append(getattr(ActsPythonBindings._propagator, f"{prefix}Propagator"))
    _propagation_algorithms.append(
        getattr(ActsPythonBindings._examples, f"{prefix}PropagationAlgorithm")
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
    def __call__(level, *args, **kwargs):
        config_or_propagator = kwargs.get("propagator", kwargs.get("config"))
        if config_or_propagator is None:
            if len(args) < 1:
                raise TypeError("Missing propagator or config argument")
            config_or_propagator = args[0]
            args = args[1:]
        for prop, alg in zip(_propagators, _propagation_algorithms):
            if isinstance(config_or_propagator, alg.Config):
                kwargs.pop("config")
                return alg(config_or_propagator, *args, **kwargs)
            if isinstance(config_or_propagator, prop):
                cfg = alg.Config(config_or_propagator)
                kwargs.pop("propagator")
                _kwargs = {}
                for k, v in kwargs.items():
                    if hasattr(cfg, k):
                        setattr(cfg, k, v)
                    else:
                        _kwargs[k] = v
                return alg(cfg, level, *args, **_kwargs)

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


_patch_config(ActsPythonBindings._examples)
