#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run a test calculation on localhost.

Usage: ./solvate.py
"""
from os import path
import click
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
from aiida_gromacs import helpers


def test_run(gromacs_code):
    """Run solvate calculation on the localhost computer.

    Uses test helpers to create AiiDA Code on the fly.
    """
    if not gromacs_code:
        # get code
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point='gromacs',
                                        computer=computer)

    # Prepare input parameters
    SolvateParameters = DataFactory('gromacs.solvate')
    parameters = SolvateParameters({'cs': 'spc216.gro',
                                    'o': '1AKI_solvated.gro'
                                    })

    SinglefileData = DataFactory('singlefile')
    grofile = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), '1AKI_newbox.gro'))
    topfile = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), '1AKI_topology.top'))

    # set up calculation
    inputs = {
        'code': gromacs_code,
        'parameters': parameters,
        'grofile': grofile,
        'topfile': topfile,
        'metadata': {
            'description': 'solvate job submission with the aiida_gromacs plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    # future = submit(CalculationFactory('gromacs'), **inputs)
    result = engine.run(CalculationFactory('gromacs.solvate'), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code):
    """Run example.

    Example usage: $ ./solvate.py --code gmx@localhost

    Alternative (creates gmx@localhost code): $ ./solvate.py

    Help: $ ./solvate.py --help
    """
    test_run(code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
