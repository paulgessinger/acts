import os

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
