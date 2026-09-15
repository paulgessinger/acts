atlas_subdir(AmvfExport)
atlas_add_component(AmvfExport src/*.cxx LINK_LIBRARIES
    AthenaBaseComps StoreGateLib BeamSpotConditionsData InDetTrackSelectionToolLib
    xAODEventInfo xAODTracking TrkParameters TrkSurfaces)
atlas_install_python_modules(python/*.py)
