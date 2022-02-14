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
    ionsmdp = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), 'ions.mdp'))
    minmdp = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), 'min.mdp'))
    eqnvt = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), 'nvt.mdp'))
    eqnpt = SinglefileData(file=path.join(path.dirname(path.realpath(__file__)), 'npt.mdp'))

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

    GromppParameters = DataFactory('gromacs.grompp')
    gromppionsparameters = GromppParameters({'o': '1AKI_ions.tpr'
                                            })

    GenionParameters = DataFactory('gromacs.genion')
    genionparameters = GenionParameters({'o': '1AKI_solvated_ions.gro',
                                         'pname': 'NA',
                                         'nname': 'CL',
                                         'neutral': 'true',
                                        })

    gromppminparameters = GromppParameters({'o': '1AKI_min.tpr'
                                            })

    MdrunParameters = DataFactory('gromacs.mdrun')
    minimiseparameters = MdrunParameters({'c': '1AKI_minimised.gro',
                                          'e': '1AKI_minimised.edr',
                                          'g': '1AKI_minimised.log',
                                          'o': '1AKI_minimised.trr',
                                          'v': 'true'
                                         })

    gromppnvtparameters = GromppParameters({'o': '1AKI_nvt.tpr',
                                            'r': '1AKI_minimised.gro'
                                            })

    nvtparameters = MdrunParameters({'c': '1AKI_nvt.gro',
                                     'e': '1AKI_nvt.edr',
                                     'g': '1AKI_nvt.log',
                                     'o': '1AKI_nvt.trr',
                                     'cpo': '1AKI_nvt.cpt',
                                     'v': 'true'
                                    })

    gromppnptparameters = GromppParameters({'o': '1AKI_npt.tpr',
                                            'r': '1AKI_nvt.gro'
                                            })

    nptparameters = MdrunParameters({'c': '1AKI_npt.gro',
                                     'e': '1AKI_npt.edr',
                                     'g': '1AKI_npt.log',
                                     'o': '1AKI_npt.trr',
                                     'cpo': '1AKI_npt.cpt',
                                     'v': 'true'
                                    })

    # set up calculation
    inputs = {
        'code': gromacs_code,
        'pdb2gmxparameters': pdb2gmxparameters,
        'editconfparameters': editconfparameters,
        'solvateparameters': solvateparameters,
        'gromppionsparameters': gromppionsparameters,
        'genionparameters': genionparameters,
        'gromppminparameters': gromppminparameters,
        'minimiseparameters': minimiseparameters,
        'gromppnvtparameters': gromppnvtparameters,
        'nvtparameters': nvtparameters,
        'gromppnptparameters': gromppnptparameters,
        'nptparameters': nptparameters,
        'pdbfile': pdbfile,
        'ionsmdp': ionsmdp,
        'minmdp': minmdp,
        'nvtmdp': eqnvt,
        'nptmdp': eqnpt,
        'metadata': {
            'description': 'setup md calculation with aiida_gromacs plugin',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    #future = engine.submit(CalculationFactory('gromacs.setup'), **inputs)
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
