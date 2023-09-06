"""Various functions for searching through previous AiiDA processes 
and appending nodes from previous processes to current process nodes.
"""

import os
import shutil
import tempfile
import re
import time
import sys

from aiida import orm
from aiida.engine import run
from aiida.plugins import CalculationFactory
from aiida.orm.nodes.process.process import ProcessState

from tests import TEST_DIR

def check_output_path(output_dir):
    """Delete existing output_dir and create a new one
    
    :param output_dir: directory path to check existence of
    :type output_dir: str
    """
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        os.makedirs(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


def run_genericMD_pdb2gmx(gromacs_code):
    """Run an instance of genericMD and return the results of the calculation.
    Used as an example calculation for tests.

    :param gromacs_code: The query entries of previous processes in the AiiDA database
    :type gromacs_code: :py:class:`aiida.orm.nodes.data.code.installed.InstalledCode`
    :returns: Results from genericMD calculation
    :rtype: dict
    """

    # input files used in pdb2gmx command
    inputs = ["pdb2gmx_1AKI_clean.pdb"]
    input_files = {}
    for filename in list(inputs):
        file_path = os.path.join(TEST_DIR, "input_files", filename)
        input_files["pdbfile"] = orm.SinglefileData(file=file_path)

    # output files produced from pdb2gmx command
    output_files = [
        "pdb2gmx_1AKI_restraints.itp",
        "pdb2gmx_1AKI_topology.top",
        "pdb2gmx_1AKI_forcefield.gro",
    ]

    # full pdb2gmx command to run
    command = (
        "pdb2gmx -i pdb2gmx_1AKI_restraints.itp "
        "-o pdb2gmx_1AKI_forcefield.gro -p pdb2gmx_1AKI_topology.top "
        "-ff oplsaa -water spce -f pdb2gmx_1AKI_clean.pdb"
    )

    # set path to temp dir
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = os.path.join(TEST_DIR, temp_dir)

    # check output_dir doesn't already exist
    check_output_path(output_dir)

    # create input dictionary for calculation.
    process_inputs = {
        "code": gromacs_code,
        "command": orm.Str(command),
        "input_files": input_files,
        "output_files": orm.List(output_files),
        "metadata": {
            "label": "generic-execute",
            "description": "Run CLI job and save input and output file provenance.",
            "options": {
                "output_filename": "file.out",
                "output_dir": output_dir,
                "parser_name": "genericMD",
            },
        },
    }

    # run calculation via aiida
    result = run(CalculationFactory("genericMD"), **process_inputs)

    return result, output_dir


def format_link_label(filename: str) -> str:
    """
    From https://github.com/sphuber/aiida-shell/blob/master/src/aiida_shell/parsers/shell.py
    Format the link label from a given filename.
    Valid link labels can only contain alphanumeric characters and
    underscores, without consecutive underscores. So all characters
    that are not alphanumeric or an underscore are converted to
    underscores, where consecutive underscores are merged into one.

    :param filename: The filename to convert to a label.
    :returns: The link label used in AiiDA provenance graphs.
    """
    alphanumeric = re.sub("[^0-9a-zA-Z_]+", "_", filename)
    link_label = re.sub("_[_]+", "_", alphanumeric)

    return link_label


def build_query():
    """
    Uses AiiDA querybuilder to find previously run processes

    :returns: the query entries
    :rtype: :py:class:`aiida.orm.querybuilder.QueryBuilder`
    """
    qb = orm.QueryBuilder()
    qb.append(orm.ProcessNode, tag='process')
    qb.order_by({orm.ProcessNode: {"ctime": "desc"}})
    return qb


def get_prev_inputs(inputs, input_labels):
    """Checks if input labels are output labels in previous processes
    and if so, adds most recent previous nodes to the new process inputs

    :param inputs: all inputs for the current process
    :type inputs: dict
    :param input_labels: input labels of the current process to search for in
        previous processes.
    :returns: updated inputs with links to previous nodes if applicable
    :rtype: dict
    """
    qb = build_query()
    added_files = []
    if qb.count() > 0:
        for entry in qb.all():
            previous_calculation = entry[0]
            for label in previous_calculation.outputs:
                if label in input_labels and label not in added_files:
                    added_files.append(label)
                    previous_output_node = \
                        previous_calculation.outputs[f"{label}"]
                    inputs[label] = previous_output_node
    return inputs


def check_prev_process(qb):
    """Wait for a previous process to finish if it is still running. The
    process state is checked every 10 seconds for up to 5 minutes and stops 
    when the process state is set to finished. If a previous processes takes 
    longer than 5 minutes, the latest submitted process is exited.

    :param qb: The queries of previous processes in the AiiDA database
    :type qb: :py:class:`aiida.orm.querybuilder.QueryBuilder`
    """
    if qb.count() > 0:
        # Get the most recently process that was already submitted to the
        # daemon and check if it has finished, wait 10s if not.
        prev_calc = qb.first()[0]
        if prev_calc.process_state != ProcessState.EXCEPTED:  
            timeout = time.time() + 60*5 # 5 minutes from now
            while prev_calc.process_state != ProcessState.FINISHED:
                print(f"Previous process status: {prev_calc.process_state}")
                print("Waiting for previous process to finish...")
                time.sleep(10)
                if time.time() > timeout:
                    sys.exit("Wait time exceeded for previous "
                             "process to complete")
        else:
            sys.exit("Previous process did not complete successfully, "
                     "please check")



def append_prev_nodes(qb, inputs, process_inputs, INPUT_DIR):
    """Checks if previous processes exists for genericMD calcs and links the 
    most recent SinglefileData type output nodes from previous processs as 
    inputs to the new process if the file names match.

    :param qb: The query entries of previous processes in the AiiDA database
    :type qb: :py:class:`aiida.orm.querybuilder.QueryBuilder`
    :param inputs: Input files for the command to be run via AiiDA
    :type inputs: list
    :param process_inputs: All inputs for the current process to be submitted
    :type process_inputs: dict
    :param INPUT_DIR: base directory where outputted files are stored.
    :returns: Updated inputs for the current process
    :rtype: dict
    """

    # if previous processes exist then check if input files are stored as
    # previous nodes and use these nodes as inputs for new process.
    if qb.count() > 0:
        # get just the name of the file from the filepath.
        stripped_inputs = []
        for inp in inputs:  # strip input file names of any paths.
            stripped_inputs.append(inp.split("/")[-1])

        prev = {}
        prev_files = []  # list of previous files already saved.
        wait_for = []
        for entry in qb.all():
            # A previous calculation exists - use its output as input for the
            # current calculation.
            previous_calculation = entry[0]
            wait_for.append(previous_calculation)
            for label in previous_calculation.outputs:
                # Get the outputs from a previous process.
                previous_output_node = previous_calculation.outputs[f"{label}"]
                # (below 2 lines does the same as above)
                # previous_output_node = orm.load_node(
                #         previous_calculation.outputs[f"{label}"].pk)

                # check if the output node is a file.
                if isinstance(previous_output_node, orm.SinglefileData):
                    prev_output_filename = previous_output_node.get_attribute(
                        "filename"
                    ) # get filename of the node
                    # check if output file is an input for new process and
                    # hasn't already been included as an input.
                    if (
                        prev_output_filename in stripped_inputs
                        and prev_output_filename not in prev_files
                    ):
                        prev_files.append(prev_output_filename)
                        prev[
                            format_link_label(prev_output_filename)
                        ] = previous_output_node

        # save input files not found in previous nodes too.
        for filename in list(inputs):
            stripped_input = filename.split("/")[-1]
            if stripped_input not in prev_files:
                prev[format_link_label(stripped_input)] = orm.SinglefileData(
                    file=os.path.join(INPUT_DIR, filename)
                )

        # update the calculation inputs dict with new dictionary of
        # input files including nodes from previous processes.
        process_inputs["input_files"] = prev
    return process_inputs