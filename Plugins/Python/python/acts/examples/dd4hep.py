import acts

import multiprocessing
import sys
import os


# Cannot conveniently catch linker errors, so we launch a suprocess to
# try importing and see if it works in order to provide a useful error message
def _import_test():
    from acts import ActsPythonBindingsDD4hep


p = multiprocessing.Process(target=_import_test)
p.start()
p.join()
if p.exitcode != 0:
    raise RuntimeError(
        "Error encountered importing DD4hep. Likely you need to set LD_LIBRARY_PATH."
    )

from acts import ActsPythonBindingsDD4hep as _dd4hep
from acts.ActsPythonBindingsDD4hep import *

import inspect

# for m in (_dd4hep,):
#     for name, cls in inspect.getmembers(m, inspect.isclass):
#         # print(name)
#         if not hasattr(cls, "Config"):
#             continue
#         if name.endswith("Detector"):
#             cls.create = _detector_create(cls)
#         else:
#             cls.__init__ = acts._make_config_adapter(cls.__init__)
