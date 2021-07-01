from acts.ActsPythonBindings._examples import *
from acts import ActsPythonBindings
import acts
from acts._adapter import _patch_config

import inspect

_propagators = []
_concrete_propagators = []
for prefix in ("Eigen", "Atlas", "StraightLine"):
    _propagators.append(getattr(ActsPythonBindings._propagator, f"{prefix}Propagator"))
    _concrete_propagators.append(
        getattr(ActsPythonBindings._propagator, f"{prefix}ConcretePropagator")
    )


def ConcretePropagator(propagator):
    for prop, prop_if in zip(_propagators, _concrete_propagators):
        if isinstance(propagator, prop):
            return prop_if(propagator)

    raise TypeError(f"Unknown propagator {type(propagator).__name__}")


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
