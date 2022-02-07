#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run a test calculation on localhost.

Usage: ./pdb2gmx.py
"""
from os import path
import click
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
from aiida_gromacs import helpers


def test_run(gromacs_code):
    """Run pdb2gmx calculation on the localhost computer.

    Uses test helpers to create AiiDA Code on the fly.
    """
    if not gromacs_code:
        # get code
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point='gromacs',
                                        computer=computer)

    # Prepare input parameters
    Pdb2gmxParameters = DataFactory('gromacs.pdb2gmx')
    parameters = Pdb2gmxParameters({'ff': 'oplsaa',
                                    'water': 'spce',
                                    'o': '1AKI_forcfield.gro',
                                    'p': '1AKI_topology.top',
                                    'i': '1AKI_restraints.itp'
                                    })

    SinglefileData = DataFactory('singlefile')
    pdbfile = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), '1AKI_clean.pdb'))

    # set up calculation
    inputs = {
        'code': gromacs_code,
        'parameters': parameters,
        'pdbfile': pdbfile,
        'metadata': {
            'description': 'pdb2gmx job submission with the aiida_gromacs plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    # future = submit(CalculationFactory('gromacs'), **inputs)
    result = engine.run(CalculationFactory('gromacs.pdb2gmx'), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code):
    """Run example.

    Example usage: $ ./pdb2gmx.py --code gmx@localhost

    Alternative (creates gmx@localhost code): $ ./pdb2gmx.py

    Help: $ ./pdb2gmx.py --help
    """
    test_run(code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
