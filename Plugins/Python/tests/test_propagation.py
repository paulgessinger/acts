import pytest

import acts
import acts.examples


def test_navigator():
    nav = acts.Navigator()
    nav = acts.Navigator(trackingGeometry=None)


def test_steppers():
    with pytest.raises(TypeError):
        acts.examples.PropagationAlgorithm()
    with pytest.raises(TypeError):
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

        prop = acts.Propagator(stepper=s, navigator=nav)

        acts.examples.PropagationAlgorithm(level=acts.logging.INFO, propagator=prop)

    assert acts.StraightLineStepper()
