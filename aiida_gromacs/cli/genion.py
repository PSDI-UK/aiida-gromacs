#!/usr/bin/env python
"""CLI utility to run gmx genion with AiiDA.

Usage: gmx_genion --help
"""

import os

import click

from aiida import cmdline, engine
from aiida.plugins import CalculationFactory, DataFactory

from aiida_gromacs import helpers


def launch(params):
    """Run genion.

    Uses helpers to add gromacs on localhost to AiiDA on the fly.
    """

    # If code is not initialised, then setup.
    gromacs_code = params.pop("code")
    if not gromacs_code:
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point="bash", computer=computer)

    # Prepare input parameters in AiiDA formats.
    SinglefileData = DataFactory("core.singlefile")
    tprfile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("s")))
    topfile = SinglefileData(file=os.path.join(os.getcwd(), params.pop("p")))

    GenionParameters = DataFactory("gromacs.genion")
    parameters = GenionParameters(params)

    # set up calculation.
    inputs = {
        "code": gromacs_code,
        "parameters": parameters,
        "tprfile": tprfile,
        "topfile": topfile,
        "metadata": {
            "description": "genion job submission with the aiida_gromacs plugin",
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # pylint: disable=unused-variable
    future = engine.submit(CalculationFactory("gromacs.genion"), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
# Input file options
@click.option("-s", default="topol.tpr", type=str, help="Input structure file")
@click.option("-n", type=str, help="Index file")
@click.option("-p", default="topol.top", type=str, help="Topology file")
# Output file options
@click.option("-o", default="out.gro", type=str, help="Output structure file")
# Other parameter options
@click.option("-np", type=str, help="Number of positive ions")
@click.option("-pname", default="NA", type=str, help="Name of positive ion")
@click.option("-pq", type=str, help="Charge of the positive ion")
@click.option("-nn", type=str, help="Number of negative ions")
@click.option("-nname", default="CL", type=str, help="Name of negative ion")
@click.option("-nq", type=str, help="Charge of the negative ion")
@click.option("-rmin", type=str, help="Minimum distance between ions and non-solvent")
@click.option("-seed", type=str, help="Seed for random number generator (0 means generate)")
@click.option("-conc", type=str, help="Specify salt concentration (mol/liter). This will add sufficient ions to reach up to the specified concentration as computed from the volume of the cell in the input .tpr file. Overrides the -np and -nn options.")
@click.option(
    "-neutral", default="false", type=str, help="Neutralise the system with ions"
)
def cli(*args, **kwargs):
    # pylint: disable=unused-argument
    # pylint: disable=line-too-long
    """Run example.

    Example usage:

    $ gmx_genion --code gmx@localhost -s 1AKI_ions.tpr -p 1AKI_topology.top -pname NA -nname CL -neutral true -o 1AKI_solvated_ions.gro

    Alternative (automatically tried to create gmx@localhost code, but requires
    gromacs to be installed and available in your environment path):

    $ gmx_genion -s 1AKI_ions.tpr -p 1AKI_topology.top -pname NA -nname CL -neutral true -o 1AKI_solvated_ions.gro

    Help: $ gmx_genion --help
    """

    launch(kwargs)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
