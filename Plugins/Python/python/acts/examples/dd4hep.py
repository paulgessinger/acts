import multiprocessing


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

from acts._adapter import _patch_config
from acts import ActsPythonBindingsDD4hep

_patch_config(ActsPythonBindingsDD4hep)

from acts.ActsPythonBindingsDD4hep import *
