# -*- coding: utf-8 -*-
"""
Calculations provided by aiida_gromacs.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida import orm
from aiida.plugins import DataFactory

Pdb2gmxParameters = DataFactory('gromacs.pdb2gmx')


class Pdb2gmxCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the 'gmx pdb2gmx' executable.

    AiiDA plugin wrapper for converting PDB files to GRO files.
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
        spec.inputs['metadata']['options']['parser_name'].default = 'gromacs.pdb2gmx'
        spec.input('metadata.options.output_filename', valid_type=str, default='aiida.out')
        spec.input('pdbfile', valid_type=orm.SinglefileData, help='Input structure (pdb).')
        spec.input('parameters', valid_type=Pdb2gmxParameters, help='Command line parameters for gmx pdb2gmx')

        spec.input('outputfile', valid_type=orm.Str, required=False, default=lambda: orm.Str('conf.gro'), help='Output forcefield compliant file (conf.gro)')
        spec.input('topfile', valid_type=orm.Str, required=False, default=lambda: orm.Str('topol.top'), help='Output (topol.top)')
        spec.input('itpfile', valid_type=orm.Str, required=False, default=lambda: orm.Str('posre.itp'), help='Output (posre.itp)')

        spec.output('stdout', valid_type=orm.SinglefileData, help='stdout')
        spec.output('outputfile', valid_type=orm.SinglefileData, help='Output forcefield compliant file (conf.gro)')
        spec.output('topfile', valid_type=orm.SinglefileData, help='Output forcefield compliant file (conf.gro)')
        spec.output('itpfile', valid_type=orm.SinglefileData, help='Output forcefield compliant file (conf.gro)')

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
            pdbfile=self.inputs.pdbfile.filename, 
            outputfile=self.inputs.outputfile.value)
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.withmpi = self.inputs.metadata.options.withmpi

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = [
            (self.inputs.pdbfile.uuid, self.inputs.pdbfile.filename, self.inputs.pdbfile.filename),
        ]
        calcinfo.retrieve_list = [self.metadata.options.output_filename,
                                  self.inputs.outputfile.value,
                                  self.inputs.topfile.value,
                                  self.inputs.itpfile.value]

        return calcinfo
