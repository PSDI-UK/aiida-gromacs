""" Tests for calculations

"""
import os

from aiida.engine import run
from aiida.plugins import CalculationFactory, DataFactory

from . import TEST_DIR


def test_process(gromacs_code):
    """Test running a grompp calculation.
    Note: this does not test that the expected outputs are created of output parsing"""

    # Prepare input parameters
    GromppParameters = DataFactory("gromacs.grompp")
    parameters = GromppParameters({"o": "grompp_1AKI_ions.tpr"})

    SinglefileData = DataFactory("core.singlefile")
    mdpfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "grompp_ions.mdp")
    )
    grofile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "grompp_1AKI_solvated.gro")
    )
    topfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "grompp_1AKI_topology.top")
    )

    # set up calculation
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "mdpfile": mdpfile,
        "grofile": grofile,
        "topfile": topfile,
        "metadata": {
            "description": "grompp test",
        },
    }

    result = run(CalculationFactory("gromacs.grompp"), **inputs)

    assert "stdout" in result
    assert "tprfile" in result
