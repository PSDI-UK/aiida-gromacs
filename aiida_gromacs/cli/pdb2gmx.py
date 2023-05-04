#!/usr/bin/env python
"""CLI utility to run gmx pdb2gmx with AiiDA.

Usage: gmx_pdb2gmx --help
"""

import click
import os
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
from aiida_gromacs import helpers


def launch(gromacs_code, pdbfile, params):
    """Run pdb2gmx.

    Uses helpers to add gromacs on localhost to AiiDA on the fly.
    """
    if not gromacs_code:
        # get code
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point='gromacs',
                                        computer=computer)

    # Prepare input parameters
    Pdb2gmxParameters = DataFactory('gromacs.pdb2gmx')
    parameters = Pdb2gmxParameters(params)

    SinglefileData = DataFactory('core.singlefile')
    pdbfile = SinglefileData(file=os.path.join(os.getcwd(), pdbfile))

    # set up calculation
    inputs = {
        'code': gromacs_code,
        'parameters': parameters,
        'pdbfile': pdbfile,
        'metadata': {
            'description': 'record pdb2gmx data provenance via the aiida_gromacs plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    future = engine.submit(CalculationFactory('gromacs.pdb2gmx'), **inputs)
    #result = engine.run(CalculationFactory('gromacs.pdb2gmx'), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
@click.option('-f', default='protein.pdb', type=str, help="Input structure file")
@click.option('-ff', required=True, type=str, help="Forcefield")
@click.option('-water', required=True, type=str, help="Water model")
@click.option('-o', default='conf.gro', type=str, help="Output structure file")
@click.option('-p', default='topol.gro', type=str, help="Topology file")
@click.option('-i', default='posre.itp', type=str, help="Include file for topology")
def cli(code, f, ff, water, o, p, i):
    """Run example.

    Example usage: 
    
    $ gmx_pdb2gmx --code gmx@localhost -f 1AKI_clean.pdb -ff oplsaa -water spce -o 1AKI_forcefield.gro -p 1AKI_topology.top -i 1AKI_restraints.itp

    Alternative (automatically tried to create gmx@localhost code, but requires
    gromacs to be installed and available in your environment path): 
    
    $ gmx_pdb2gmx -f 1AKI_clean.pdb -ff oplsaa -water spce -o 1AKI_forcefield.gro -p 1AKI_topology.top -i 1AKI_restraints.itp

    Help: $ gmx_pdb2gmx --help
    """
    
    # Place CLI params in a dict.
    params={'ff': ff,
            'water': water,
            'o': o,
            'p': p,
            'i': i
           }
    
    launch(code, f, params)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter


