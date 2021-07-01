import pytest

import acts
import acts.examples


def test_navigator(conf_const):
    nav = conf_const(acts.Navigator)
    nav = conf_const(acts.Navigator, trackingGeometry=None)


def test_steppers(conf_const):
    with pytest.raises(TypeError):
        acts.examples.PropagationAlgorithm()
    with pytest.raises(ValueError):
        acts.examples.PropagationAlgorithm(level=acts.logging.INFO)

    with pytest.raises(TypeError):
        acts.Propagator()

    nav = acts.Navigator()

    with pytest.raises(TypeError):
        acts.Propagator(navigator=nav)

    for stepper in (acts.EigenStepper, acts.AtlasStepper):
        with pytest.raises(TypeError):
            stepper()
        s = stepper(acts.NullBField())
        assert s

        prop = acts.examples.ConcretePropagator(
            acts.Propagator(stepper=s, navigator=nav)
        )

        # cfg = acts.examples.PropagationAlgorithm.Config()
        # cfg.propagatorImpl = prop
        # acts.examples.PropagationAlgorithm(level=acts.logging.INFO, config=cfg)
        # acts.examples.PropagationAlgorithm(cfg, acts.logging.INFO)

        assert conf_const(
            acts.examples.PropagationAlgorithm,
            level=acts.logging.INFO,
            propagatorImpl=prop,
        )

    assert acts.StraightLineStepper()
