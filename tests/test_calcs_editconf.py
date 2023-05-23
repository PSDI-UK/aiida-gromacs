""" Tests for calculations

"""
import os

from aiida.engine import run
from aiida.plugins import CalculationFactory, DataFactory

from . import TEST_DIR


def test_process(gromacs_code):
    """Test running a editconf calculation.
    Note: this does not test that the expected outputs are created of output parsing"""

    # Prepare input parameters
    EditconfParameters = DataFactory("gromacs.editconf")
    parameters = EditconfParameters(
        {"center": "0", "d": "1.0", "bt": "cubic", "o": "editconf_1AKI_newbox.gro"}
    )

    SinglefileData = DataFactory("core.singlefile")
    grofile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "editconf_1AKI_forcefield.gro")
    )

    # set up calculation
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "grofile": grofile,
        "metadata": {
            "description": "editconf test",
        },
    }

    result = run(CalculationFactory("gromacs.editconf"), **inputs)

    assert "stdout" in result
    assert "grofile" in result
