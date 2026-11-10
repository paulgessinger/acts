import acts
import pytest

from helpers import detrayEnabled

u = acts.UnitConstants


@pytest.mark.skipif(not detrayEnabled, reason="Detray plugin not available")
def test_detray_itk_bindings():
    from acts import detray
    import acts.examples.detray as exdetray

    # The ITk detector type sits next to the ODD one
    assert hasattr(detray, "DetrayDetectorITk")
    assert hasattr(detray, "readITk")
    assert hasattr(detray.DetrayGeometryConverter, "convertITk")
    for tag in ("ODD", "ITk"):
        assert hasattr(exdetray, f"StraightLinePropagator{tag}")
        assert hasattr(exdetray, f"HelixPropagator{tag}")

    if hasattr(exdetray, "NavigationValidationConfig"):
        cfg = exdetray.NavigationValidationConfig()
        assert cfg.rays and cfg.helices
        assert not cfg.displayOnFailure
        cfg.nTracks = 10
        cfg.pTRange = [1 * u.GeV, 2 * u.GeV]
        assert cfg.nTracks == 10


@pytest.mark.skipif(not detrayEnabled, reason="Detray plugin not available")
def test_detray_surface_grid_policy():
    from acts import detray
    import acts.json

    cfg = detray.DetraySurfaceGridConfig()
    assert cfg.vertexZHalfRange == pytest.approx(150 * u.mm)
    assert cfg.minPt == pytest.approx(0.5 * u.GeV)
    assert cfg.numberOfBinsFactor == pytest.approx(2.0)
    assert cfg.lookupOnLayer
    cfg.minPt = 1 * u.GeV
    cfg.numberOfBinsFactor = None
    assert cfg.numberOfBinsFactor is None

    factory = detray.addDetraySurfaceGridPolicy(
        acts.NavigationPolicyFactory.make().add(acts.TryAllNavigationPolicy), cfg
    )
    assert factory is not None

    payloadCfg = detray.DetrayPayloadConverter.Config()
    assert payloadCfg.surfaceGridFill is not None
    payloadCfg.surfaceGridFill = None
    assert payloadCfg.surfaceGridFill is None
    payloadCfg.surfaceGridFill = cfg
    assert payloadCfg.surfaceGridFill.minPt == pytest.approx(1 * u.GeV)

    jsonCfg = acts.json.TrackingGeometryJsonConverter.Config()
    detray.addDetrayNavigationPolicies(jsonCfg)
    acts.json.TrackingGeometryJsonConverter(jsonCfg)
