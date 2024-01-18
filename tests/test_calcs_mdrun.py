""" Tests for calculations

"""
import os

from aiida.engine import run
from aiida.orm import Dict
from aiida.plugins import CalculationFactory, DataFactory

from . import TEST_DIR


def run_mdrun(gromacs_code):
    """Run an instance of mdrun and return the results."""

    # Prepare input parameters
    MdrunParameters = DataFactory("gromacs.mdrun")
    parameters = MdrunParameters(
        {
            "c": "mdrun_1AKI_minimised.gro",
            "e": "mdrun_1AKI_minimised.edr",
            "g": "mdrun_1AKI_minimised.log",
            "o": "mdrun_1AKI_minimised.trr",
            "v": "true",
            "ntomp": "5",
            "ntmpi": "1",
        }
    )

    SinglefileData = DataFactory("core.singlefile")
    tprfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "mdrun_1AKI_em.tpr")
    )

    # set up calculation
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "tprfile": tprfile,
        "metadata": {
            "description": "mdrun test",
        },
    }

    result = run(CalculationFactory("gromacs.mdrun"), **inputs)

    return result


def test_process(gromacs_code):
    """Test running a mdrun calculation.
    Note: this does not test that the expected outputs are created of output parsing"""

    result = run_mdrun(gromacs_code)

    assert "stdout" in result
    assert "trrfile" in result
    assert "grofile" in result
    assert "logfile" in result
    assert "enfile" in result
    assert "logfile_metadata" in result


def test_file_name_match(gromacs_code):
    """Test that the file names returned match what was specified on inputs."""

    result = run_mdrun(gromacs_code)

    assert result["stdout"].base.repository.list_object_names()[0] == "mdrun.out"
    assert (
        result["trrfile"].base.repository.list_object_names()[0]
        == "mdrun_1AKI_minimised.trr"
    )
    assert (
        result["grofile"].base.repository.list_object_names()[0]
        == "mdrun_1AKI_minimised.gro"
    )
    assert (
        result["logfile"].base.repository.list_object_names()[0]
        == "mdrun_1AKI_minimised.log"
    )
    assert (
        result["enfile"].base.repository.list_object_names()[0]
        == "mdrun_1AKI_minimised.edr"
    )
    assert isinstance(result["logfile_metadata"], Dict)
