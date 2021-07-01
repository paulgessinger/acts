import pytest

import acts


def test_material(conf_const):
    with pytest.raises(TypeError):
        acts.examples.RootMaterialDecorator()
    fileName = "blubb.root"
    try:
        conf_const(
            acts.examples.RootMaterialDecorator,
            level=acts.logging.INFO,
            fileName=fileName,
        )
    except RuntimeError as e:
        assert fileName in str(e)
