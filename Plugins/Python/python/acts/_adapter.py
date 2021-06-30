import inspect
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


def _make_config_constructor(cls):
    fn = cls.__init__

    @functools.wraps(fn)
    def wrapped(self, *args, **kwargs):
        _kwargs = {}
        for k in list(kwargs.keys()):
            if hasattr(cls, k):
                _kwargs[k] = kwargs.pop(k)

        fn(self, *args, **kwargs)

        for k, v in _kwargs.items():
            setattr(self, k, v)

    return wrapped


def _patch_config_constructor(cls):
    cls.__init__ = _make_config_constructor(cls)


def _patch_config(m):
    for name, cls in inspect.getmembers(m, inspect.isclass):
        if hasattr(cls, "Config"):
            cls.__init__ = _make_config_adapter(cls.__init__)

        if name == "Config":
            _patch_config_constructor(cls)
