#!/usr/bin/env python
"""
Launch a calculation using the 'general-MD'
"""

import os
from pathlib import Path
import re
from time import sleep

import click

from aiida import cmdline, engine, orm
from aiida.common import exceptions
from aiida.orm.nodes.process.process import ProcessState
from aiida.orm.querybuilder import QueryBuilder
from aiida.plugins import CalculationFactory
from aiida import load_profile

# from aiida.common.exceptions import NotExistent
from aiida_gromacs import helpers

# set base path for input files.
INPUT_DIR = os.path.join(os.getcwd())
profile = load_profile()
computer = helpers.get_computer()
code = helpers.get_code(entry_point="gromacs", computer=computer)
code = helpers.get_code(entry_point="bash", computer=computer)

def format_link_label(filename: str) -> str:
    """
    From: https://github.com/sphuber/aiida-shell/blob/master/src/aiida_shell/parsers/shell.py
    Format the link label from a given filename.
    Valid link labels can only contain alphanumeric characters and
        underscores, without consecutive underscores. So all characters
        that are not alphanumeric or an underscore are converted to
        underscores, where consecutive underscores are merged into one.
    Additional: Label cannot start with a number or underscore.
    :param filename: The filename.
    :returns: The link label.
    """
    alphanumeric = re.sub("[^0-9a-zA-Z_]+", "_", filename)
    link_label = re.sub("_[_]+", "_", alphanumeric)

    return link_label


def check_prev_process(qb):
    """Wait for previous process to finish if running"""
    if qb.count() > 0:
        # Get the most recently process that was already submitted to the
        # daemon and check if it has finished, wait 10s if not.
        prev_calc = qb.first()[0]

        while prev_calc.process_state != ProcessState.FINISHED:
            print(f"Previous process status: {prev_calc.process_state}")
            print("Waiting for previous process to finish...")
            sleep(10)


def append_prev_nodes(qb, inputs, inputs_):
    """If previous processes exist, link nodes to files to new process"""

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
            # Get the outputs from a previous process.
            # previous_output = previous_calculation.outputs
            # previous_node = orm.load_node(previous_calculation)

            for label in previous_calculation.outputs:
                previous_output_node = previous_calculation.outputs[f"{label}"]
                # (below 2 lines does the same as above)
                # previous_output_node = orm.load_node(
                #         previous_calculation.outputs[f"{label}"].pk)
                if isinstance(previous_output_node, orm.SinglefileData):
                    # check if the output node is a file.
                    prev_output_filename = previous_output_node.get_attribute(
                        "filename"
                    )
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
        inputs_["input_files"] = prev
    return inputs_


def launch_generalMD(options):
    """Run general-MD"""

    code = options["code"]
    command = options["command"]
    inputs = options["inputs"]
    outputs = options["outputs"]
    output_dir = options["output_dir"]

    print(f"command: {command}")

    if not code:
        raise exceptions.NonExistent("Code has not been set.")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    MyAppCalculation = CalculationFactory("general-MD")

    # Check if a previous calculation with the same input parameter
    # value has been stored by loading the QueryBuilder and append
    # all previous jobs ordered by newest first.
    qb = QueryBuilder()
    qb.append(MyAppCalculation, tag="calcjob")
    qb.order_by({MyAppCalculation: {"ctime": "desc"}})

    check_prev_process(qb)

    # Save list of input files to a dict with keys that are formatted
    # file names and values that are SinglefileData.
    input_files = {}
    for filename in list(inputs):
        file_path = os.path.join(INPUT_DIR, filename)
        stripped_input = filename.split("/")[-1]
        input_files[format_link_label(stripped_input)] = orm.SinglefileData(
            file=file_path
        )

    # Keep the output filenames as a list.
    output_files = list(outputs)

    # create input dictionary for calculation.
    inputs_ = {
        "code": code,
        "command": orm.Str(command),
        "input_files": input_files,
        "output_files": orm.List(output_files),
        "metadata": {
            "label": "general-execute",
            "description": "Run CLI job and save input and output " "file provenance.",
            "options": {
                "output_filename": "file.log",
                "output_dir": output_dir,
                "parser_name": "general-MD",
            },
        },
    }

    if qb.count() > 0:
        inputs_ = append_prev_nodes(qb, inputs, inputs_)

    # Submit your calculation to the aiida daemon
    # pylint: disable=unused-variable
    future = engine.submit(CalculationFactory("general-MD"), **inputs_)
    # future = engine.submit(process)
    print(f"Submitted calculation: {future}\n")


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
# @click.option(
#     "--code",
#     type=str,
#     multiple=True,
#     help="The code used to in the executable.",
# )
@click.option(
    "--command",
    type=str,
    help="The full command used to run the job enclosed in quotes.",
)
@click.option(
    "--inputs",
    type=str,
    multiple=True,
    help="Input file name used in the command. "
    "Include the local path to these files.",
)
@click.option(
    "--outputs", multiple=True, type=str, help="Output file name used in the command."
)
@click.option(
    "--output_dir",
    default=os.path.join(os.getcwd()) + "/outputs",
    type=str,
    help="Absolute path of directory where files are saved.",
)
def cli(**kwargs):
    """Run general-MD

    Example usage:

    $ ./launch.py --code gmx@localhost
    --command "pdb2gmx -i 1AKI_restraints.itp -o 1AKI_forcfield.gro
    -p 1AKI_topology.top -ff oplsaa -water spce -f 1AKI_clean.pdb"
    --inputs 1AKI_clean.pdb
    --outputs 1AKI_restraints.itp
    --outputs 1AKI_topology.top
    --outputs 1AKI_forcfield.gro

    Help: $ ./launch.py --help
    """
    launch_generalMD(kwargs)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
