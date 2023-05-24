""" Tests for calculations

"""
import os

from aiida.engine import run
from aiida.plugins import CalculationFactory, DataFactory

from . import TEST_DIR


def test_process(bash_code):
    """Test running a genion calculation.
    Note: this does not test that the expected outputs are created of output parsing"""

    # Prepare input parameters
    GenionParameters = DataFactory("gromacs.genion")
    parameters = GenionParameters(
        {
            "o": "genion_1AKI_solvated_ions.gro",
            "pname": "NA",
            "nname": "CL",
            "neutral": "true",
        }
    )

    SinglefileData = DataFactory("core.singlefile")
    tprfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "genion_1AKI_ions.tpr")
    )
    topfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "genion_1AKI_topology.top")
    )

    # set up calculation
    inputs = {
        "code": bash_code,
        "parameters": parameters,
        "tprfile": tprfile,
        "topfile": topfile,
        "metadata": {
            "description": "genion test",
        },
    }

    result = run(CalculationFactory("gromacs.genion"), **inputs)

    assert "stdout" in result
    assert "grofile" in result
    assert "topfile" in result
