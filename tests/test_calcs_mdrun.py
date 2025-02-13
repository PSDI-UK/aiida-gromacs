""" Tests for calculations

"""
import os

from aiida.engine import run
from aiida.orm import Dict
from aiida.plugins import CalculationFactory, DataFactory

from aiida_gromacs.data.plumed_input import PlumedInputData

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
            # "ntmpi": "1", # turn off omp and mpi for gmx patched with plumed
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


def run_mdrun_plumed(gromacs_code):
    """Run an instance of mdrun and return the results."""

    # Prepare input parameters
    MdrunParameters = DataFactory("gromacs.mdrun")
    parameters = MdrunParameters(
        {
            "c": "plumed_mdrun_prod.gro",
            "e": "plumed_mdrun_prod.edr",
            "g": "plumed_mdrun_prod.log",
            "o": "plumed_mdrun_prod.trr",
            # "plumed": "plumed_mdrun_prod.dat",
            "v": "true",
            "ntomp": "5",
            # "ntmpi": "1", # turn off omp and mpi for gmx patched with plumed
        }
    )

    SinglefileData = DataFactory("core.singlefile")
    tprfile = SinglefileData(
        file=os.path.join(TEST_DIR, "input_files", "plumed_mdrun_prod.tpr")
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

    # Set the plumed script as a PlumedInputData type node
    inputs["plumed_file"] = PlumedInputData(
        file=os.path.join(TEST_DIR, "input_files", "plumed_mdrun_prod.dat")
    )
    # Find the inputs and outputs referenced in the plumed script
    calc_inputs, calc_outputs = inputs["plumed_file"].calculation_inputs_outputs
    # add input files and dirs referenced in plumed file into inputs
    inputs.update(calc_inputs)
    inputs.update(calc_outputs)

    result = run(CalculationFactory("gromacs.mdrun"), **inputs)

    return result


def test_process_plumed(gromacs_code):
    """Test running a mdrun calculation.
    Note: this does not test that the expected outputs are created of output parsing"""

    result = run_mdrun_plumed(gromacs_code)

    assert "stdout" in result
    assert "trrfile" in result
    assert "grofile" in result
    assert "logfile" in result
    assert "enfile" in result
    assert "logfile_metadata" in result
    assert "plumed_HILLS" in result
    assert "plumed_COLVAR" in result


def test_file_name_match_plumed(gromacs_code):
    """Test that the file names returned match what was specified on inputs."""

    result = run_mdrun_plumed(gromacs_code)

    assert result["stdout"].base.repository.list_object_names()[0] == "mdrun.out"
    assert (
        result["trrfile"].base.repository.list_object_names()[0]
        == "plumed_mdrun_prod.trr"
    )
    assert (
        result["grofile"].base.repository.list_object_names()[0]
        == "plumed_mdrun_prod.gro"
    )
    assert (
        result["logfile"].base.repository.list_object_names()[0]
        == "plumed_mdrun_prod.log"
    )
    assert (
        result["enfile"].base.repository.list_object_names()[0]
        == "plumed_mdrun_prod.edr"
    )
    assert isinstance(result["logfile_metadata"], Dict)
    assert result["plumed_HILLS"].base.repository.list_object_names()[0] == "HILLS"
    assert result["plumed_COLVAR"].base.repository.list_object_names()[0] == "COLVAR"
