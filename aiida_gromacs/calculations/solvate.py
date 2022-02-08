# -*- coding: utf-8 -*-
"""
Calculations provided by aiida_gromacs.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida import orm
from aiida.plugins import DataFactory

SolvateParameters = DataFactory('gromacs.solvate')


class SolvateCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the 'gmx solvate' executable.

    AiiDA plugin wrapper for solvating a molecular system.
    """

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        # yapf: disable
        super().define(spec)

        # set default values for AiiDA options
        spec.inputs['metadata']['options']['resources'].default = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1,
        }
        spec.inputs['metadata']['options']['parser_name'].default = 'gromacs.solvate'
        spec.input('metadata.options.output_filename', valid_type=str, default='aiida.out')
        spec.input('grofile', valid_type=orm.SinglefileData, help='Input structure')
        spec.input('topfile', valid_type=orm.SinglefileData, help='Input topology')
        spec.input('parameters', valid_type=SolvateParameters, help='Command line parameters for gmx solvate.')

        spec.output('stdout', valid_type=orm.SinglefileData, help='stdout')
        spec.output('outputfile', valid_type=orm.SinglefileData, help='Output forcefield compliant file.')
        spec.output('topfile', valid_type=orm.SinglefileData, help='Output forcefield compliant file.')

        spec.exit_code(300, 'ERROR_MISSING_OUTPUT_FILES', message='Calculation did not produce all expected output files.')

    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files
            needed by the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = self.inputs.parameters.cmdline_params(
            grofile=self.inputs.grofile.filename,
            topfile=self.inputs.topfile.filename)
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.withmpi = self.inputs.metadata.options.withmpi

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = [
            (self.inputs.grofile.uuid, self.inputs.grofile.filename, self.inputs.grofile.filename), (self.inputs.topfile.uuid, self.inputs.topfile.filename, self.inputs.topfile.filename),
        ]
        calcinfo.retrieve_list = [self.metadata.options.output_filename,
                                  self.inputs.parameters['o'],
                                  self.inputs.topfile.filename]

        return calcinfo
