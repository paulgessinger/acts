# Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
"""Attach a parameter exporter to the ART reconstruction's vertex inputs."""
import os
from AthenaConfiguration.ComponentFactory import CompFactory


def addExport(flags, cfg):
    from InDetTrackSelectionTool.InDetTrackSelectionToolConfig import (
        VtxInDetTrackSelectionCfg,
    )
    from BeamSpotConditions.BeamSpotConditionsConfig import BeamSpotCondAlgCfg

    cfg.merge(BeamSpotCondAlgCfg(flags))
    selector = cfg.popToolsAndMerge(
        VtxInDetTrackSelectionCfg(flags, name="AmvfExportSelection")
    )
    cfg.addEventAlgo(
        CompFactory.AmvfTrackExport(
            "AmvfTrackExport",
            TrackSelector=selector,
            Tracks=os.environ.get("AMVF_TRACK_COLLECTION", "InDetTrackParticles"),
            OutputDirectory=os.environ.get("AMVF_EXPORT_DIR", "amvf-tracks"),
            UseBeamConstraint=flags.Tracking.PriVertex.useBeamConstraint,
        )
    )
    # Save effective flags, including vertex and selection settings.
    with open("amvf-config-flags.txt", "w") as output:
        import contextlib

        with contextlib.redirect_stdout(output):
            flags.dump(
                pattern="Tracking.*|Acts.*|GeoModel.*|IOVDb.*|Beam.*|Input.*",
                evaluate=True,
            )
    with open("amvf-reco-config.pkl", "wb") as output:
        cfg.store(output)
