from acts._acts._examples import *
import acts._acts._examples

from acts import _acts

_propagators = []
_propagation_algorithms = []
for prefix in ("Eigen", "Atlas", "StraightLine"):
    _propagators.append(getattr(acts._acts._propagator, f"{prefix}Propagator"))
    _propagation_algorithms.append(
        getattr(acts._acts._examples, f"{prefix}PropagationAlgorithm")
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


import functools


def _make_config_adapter(fn):
    @functools.wraps(fn)
    def wrapped(self, *args, **kwargs):
        if len(args) > 0:
            maybe_config = args[0]
            if isinstance(maybe_config, type(self).Config):
                # is already config, nothing to do here
                fn(self, maybe_config, *args[1:], **kwargs)
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


import inspect

for name, cls in inspect.getmembers(acts._acts._examples, inspect.isclass):
    # print(name)
    if not hasattr(cls, "Config"):
        continue
    if name.endswith("Detector"):
        cls.create = _detector_create(cls)
    else:
        cls.__init__ = _make_config_adapter(cls.__init__)
