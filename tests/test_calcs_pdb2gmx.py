""" Tests for calculations

"""
import os

from aiida.engine import run
from aiida.plugins import CalculationFactory, DataFactory

from . import TEST_DIR


def test_process(gromacs_code):
    """Test running a pdb2gmx calculation.
    Note: this does not test that the expected outputs are created of output parsing"""

    # Prepare input parameters
    DiffParameters = DataFactory("gromacs.pdb2gmx")
    parameters = DiffParameters(
        {
            "ff": "oplsaa",
            "water": "spce",
            "o": "1AKI_forcefield.gro",
            "p": "1AKI_topology.top",
            "i": "1AKI_restraints.itp",
        }
    )

    SinglefileData = DataFactory("core.singlefile")
    pdbfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "1AKI_clean.pdb")
    )

    # set up calculation
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "pdbfile": pdbfile,
        "metadata": {
            "description": "pdb2gmx test",
        },
    }

    result = run(CalculationFactory("gromacs.pdb2gmx"), **inputs)

    assert "stdout" in result
    assert "grofile" in result
    assert "itpfile" in result
    assert "topfile" in result
