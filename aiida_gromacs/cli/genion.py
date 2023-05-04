#!/usr/bin/env python
"""CLI utility to run gmx genion with AiiDA.

Usage: gmx_genion --help
"""

import click
import os
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
from aiida_gromacs import helpers


def launch(gromacs_code, tprfile, topfile, params):
    """Run genion.

    Uses helpers to add gromacs on localhost to AiiDA on the fly.
    """
    if not gromacs_code:
        # get code
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point='bash',
                                        computer=computer)

    # Prepare input parameters
    GenionParameters = DataFactory('gromacs.genion')
    parameters = GenionParameters(params)

    SinglefileData = DataFactory('core.singlefile')
    tprfile = SinglefileData(file=os.path.join(os.getcwd(), tprfile))
    topfile = SinglefileData(file=os.path.join(os.getcwd(), topfile))

    # set up calculation
    inputs = {
        'code': gromacs_code,
        'parameters': parameters,
        'tprfile': tprfile,
        'topfile': topfile,
        'metadata': {
            'description': 'genion job submission with the aiida_gromacs plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    future = submit(CalculationFactory('gromacs'), **inputs)
    #result = engine.run(CalculationFactory('gromacs.genion'), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
@click.option('-s', default='topol.tpr', type=str, help="Input structure file")
@click.option('-p', default='topol.top', type=str, help="Topology file")
@click.option('-pname', default='NA', type=str, help="Name of positive ion")
@click.option('-nname', default='CL', type=str, help="Name of negative ion")
@click.option('-neutral', default='false', type=str, help="Neutralise the system with ions")
@click.option('-o', default='out.gro', type=str, help="Output structure file")
def cli(code, s, p, pname, nname, neutral, o):
    """Run example.

    Example usage: 
    
    $ gmx_genion --code gmx@localhost -s 1AKI_ions.tpr -p 1AKI_topology.top -pname NA -nname CL -neutral true -o 1AKI_solvated_ions.gro

    Alternative (automatically tried to create gmx@localhost code, but requires
    gromacs to be installed and available in your environment path): 
    
    $ gmx_genion -s 1AKI_ions.tpr -p 1AKI_topology.top -pname NA -nname CL -neutral true -o 1AKI_solvated_ions.gro

    Help: $ gmx_genion --help
    """

    # Place CLI params in a dict.
    params={'pname': pname,
            'nname': nname,
            'neutral': neutral,
            'o': o,
           }

    launch(code, s, p, params)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
