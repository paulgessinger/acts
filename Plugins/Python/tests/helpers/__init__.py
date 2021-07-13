import os
from acts.examples import BareAlgorithm

geant4Enabled = any(v.startswith("G4") for v in os.environ.keys())
rootEnabled = "ROOTSYS" in os.environ

if rootEnabled:
    try:
        import ROOT
    except ImportError:
        import warnings

        warnings.warn(
            "ROOT likely built without/with incompatible PyROOT. Skipping tests that need ROOT"
        )
        rootEnabled = False

dd4hepEnabled = "DD4hep_DIR" in os.environ

try:
    import acts.examples.hepmc3

    hepmc3Enabled = True
except ImportError:
    hepmc3Enabled = False


class AssertCollectionExistsAlg(BareAlgorithm):
    events_seen = 0

    def __init__(self, collection, *args, **kwargs):
        self.collection = collection
        BareAlgorithm.__init__(self, *args, **kwargs)

    def execute(self, ctx):
        assert ctx.eventStore.exists(
            self.collection
        ), f"{self.collection} does not exist"
        self.events_seen += 1
        return acts.examples.ProcessCode.SUCCESS
