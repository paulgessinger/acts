import os

geant4Enabled = any(v.startswith("G4") for v in os.environ.keys())
rootEnabled = "ROOTSYS" in os.environ
dd4hepEnabled = True  # @TODO: make this a check
