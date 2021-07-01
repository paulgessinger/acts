import pytest


def kwargsConstructor(cls, *args, **kwargs):
    return cls(*args, **kwargs)


def configKwConstructor(cls, *args, **kwargs):
    assert hasattr(cls, "Config")
    _kwargs = {}
    if "level" in kwargs:
        _kwargs["level"] = kwargs.pop("level")
    config = cls.Config()
    for k, v in kwargs.items():
        setattr(config, k, v)
    return cls(*args, config=config, **_kwargs)


def configPosConstructor(cls, *args, **kwargs):
    assert hasattr(cls, "Config")
    _kwargs = {}
    if "level" in kwargs:
        _kwargs["level"] = kwargs.pop("level")
    config = cls.Config()
    for k, v in kwargs.items():
        setattr(config, k, v)

    return cls(config, *args, **_kwargs)


@pytest.fixture(params=[configPosConstructor, configKwConstructor, kwargsConstructor])
def conf_const(request):
    return request.param
