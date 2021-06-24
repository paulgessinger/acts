from ._acts import *
from ._acts import __version__
from . import _acts

import inspect
import functools


def Propagator(stepper, nav):
    for prefix in ("Eigen", "Atlas", "StraightLine"):
        _stepper = getattr(_acts, f"{prefix}Stepper")
        if isinstance(stepper, _stepper):
            return getattr(_acts._propagator, f"{prefix}Propagator")(stepper, nav)
    raise TypeError(f"Unknown stepper {type(stepper).__name__}")


def _make_config_adapter(fn):
    @functools.wraps(fn)
    def wrapped(self, *args, **kwargs):
        if len(args) > 0:
            maybe_config = args[0]
            if isinstance(maybe_config, type(self).Config):
                # is already config, nothing to do here
                fn(self, maybe_config, *args[1:], **kwargs)
                return

        if "config" in kwargs:
            config = kwargs.pop("config")
            fn(self, config, *args, **kwargs)
            return

        cfg = type(self).Config()
        _kwargs = {}
        for k, v in kwargs.items():
            if hasattr(cfg, k):
                setattr(cfg, k, v)
            else:
                _kwargs[k] = v
        fn(self, cfg, *args, **_kwargs)

    return wrapped


for m in (_acts,):
    for name, cls in inspect.getmembers(m, inspect.isclass):
        if not hasattr(cls, "Config"):
            continue

        cls.__init__ = _make_config_adapter(cls.__init__)
