#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run a test calculation on localhost.

Usage: ./editconf.py
"""
from os import path
import click
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
from aiida_gromacs import helpers


def test_run(gromacs_code):
    """Run editconf calculation on the localhost computer.

    Uses test helpers to create AiiDA Code on the fly.
    """
    if not gromacs_code:
        # get code
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point='gromacs',
                                        computer=computer)

    # Prepare input parameters
    EditconfParameters = DataFactory('gromacs.editconf')
    parameters = EditconfParameters({'center': '0',
                                    'd': '1.0',
                                    'bt': 'cubic',
                                    'o': '1AKI_newbox.gro'
                                    })

    SinglefileData = DataFactory('singlefile')
    grofile = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), '1AKI_forcefield.gro'))

    # set up calculation
    inputs = {
        'code': gromacs_code,
        'parameters': parameters,
        'grofile': grofile,
        'metadata': {
            'description': 'editconf job submission with the aiida_gromacs plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    # future = submit(CalculationFactory('gromacs'), **inputs)
    result = engine.run(CalculationFactory('gromacs.editconf'), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code):
    """Run example.

    Example usage: $ ./editconf.py --code gmx@localhost

    Alternative (creates gmx@localhost code): $ ./editconf.py

    Help: $ ./editconf.py --help
    """
    test_run(code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
