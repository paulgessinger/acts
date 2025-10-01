import logging
from rich.logging import RichHandler


def configure_logging(level: int):
    FORMAT = "%(message)s"
    logging.basicConfig(
        level=level, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
    )


class Logger:
    _logger: logging.Logger

    def __init__(self, logger: logging.Logger):
        self._logger = logger

    @staticmethod
    def patch_kwargs(kwargs: dict) -> dict:
        kwargs.setdefault("extra", {})
        kwargs["extra"].setdefault("markup", True)
        return kwargs

    def info(self, *args, **kwargs):
        self._logger.info(*args, **self.patch_kwargs(kwargs))

    def debug(self, *args, **kwargs):
        self._logger.debug(*args, **self.patch_kwargs(kwargs))

    def warning(self, *args, **kwargs):
        self._logger.warning(*args, **self.patch_kwargs(kwargs))

    def error(self, *args, **kwargs):
        self._logger.error(*args, **self.patch_kwargs(kwargs))


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)
