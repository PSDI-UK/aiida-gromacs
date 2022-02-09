#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run a test md setup workflow on localhost.

Usage: ./md-setup.py
"""
from os import path
import click
from aiida import cmdline, engine
from aiida.plugins import DataFactory, WorkflowFactory
from aiida_gromacs import helpers


def test_run(gromacs_code):
    """Run setup workflow on the localhost computer.

    Uses test helpers to create AiiDA Code on the fly.
    """
    if not gromacs_code:
        # get code
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point='gromacs',
                                        computer=computer)

    # Prepare input parameters
    SinglefileData = DataFactory('singlefile')
    pdbfile = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), '1AKI_clean.pdb'))

    Pdb2gmxParameters = DataFactory('gromacs.pdb2gmx')
    pdb2gmxparameters = Pdb2gmxParameters({'ff': 'oplsaa',
                                    'water': 'spce',
                                    'o': '1AKI_forcfield.gro',
                                    'p': '1AKI_topology.top',
                                    'i': '1AKI_restraints.itp'
                                    })

    EditconfParameters = DataFactory('gromacs.editconf')
    editconfparameters = EditconfParameters({'center': '0',
                                    'd': '1.0',
                                    'bt': 'cubic',
                                    'o': '1AKI_newbox.gro'
                                    })

    SolvateParameters = DataFactory('gromacs.solvate')
    solvateparameters = SolvateParameters({'cs': 'spc216.gro',
                                    'o': '1AKI_solvated.gro'
                                    })

    # set up calculation
    inputs = {
        'code': gromacs_code,
        'pdb2gmxparameters': pdb2gmxparameters,
        'editconfparameters': editconfparameters,
        'solvateparameters': solvateparameters,
        'pdbfile': pdbfile,
        'metadata': {
            'description': 'setup md calculation with aiida_gromacs plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    # future = submit(CalculationFactory('gromacs'), **inputs)
    result = engine.run(WorkflowFactory('gromacs.setup'), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code):
    """Run example workflow.

    Example usage: $ ./md-setup.py --code gmx@localhost

    Alternative (creates gmx@localhost code): $ ./md-setup.py

    Help: $ ./md-setup.py --help
    """
    test_run(code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
