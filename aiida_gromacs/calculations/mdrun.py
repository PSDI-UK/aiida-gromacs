"""
Calculations provided by aiida_gromacs.

This calculation configures the ability to use the 'gmx mdrun' executable.
"""
import os

from aiida.common import CalcInfo, CodeInfo
from aiida.engine import CalcJob
from aiida.orm import SinglefileData, Dict, Str
from aiida.plugins import DataFactory

MdrunParameters = DataFactory("gromacs.mdrun")


class MdrunCalculation(CalcJob):
    """
    AiiDA calculation plugin wrapping the 'gmx mdrun' executable.

    AiiDA plugin wrapper for converting PDB files to GRO files.
    """

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        # yapf: disable
        super().define(spec)

        # Define inputs and outputs of the calculation.
        spec.input('command',
                valid_type=Str, required=False,
                help='The command used to execute the job.')

        # set default values for AiiDA options
        # TODO: something changed about withmpi in aiida-2.4.0, needs investigation.
        spec.inputs['metadata']['options']['withmpi'].default = False
        # TODO: remove this for production release.
        spec.inputs['metadata']['options']['max_wallclock_seconds'].default = 86400
        spec.inputs['metadata']['options']['resources'].default = {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1,
            'num_cores_per_mpiproc': 5,
        }

        # Required inputs.
        spec.inputs['metadata']['options']['parser_name'].default = 'gromacs.mdrun'
        spec.input('metadata.options.output_filename', valid_type=str, default='mdrun.out')
        spec.input('tprfile', valid_type=SinglefileData, help='Input structure.')
        spec.input('parameters', valid_type=MdrunParameters, help='Command line parameters for gmx mdrun')
        spec.input('metadata.options.output_dir', valid_type=str, default=os.getcwd(),
                help='Directory where output files will be saved when parsed.')

        # Optional inputs.
        spec.input('cpi_file', valid_type=SinglefileData, required=False, help='Checkpoint file')
        spec.input('table_file', valid_type=SinglefileData, required=False, help='xvgr/xmgr file')
        spec.input('tableb_file', valid_type=SinglefileData, required=False, help='xvgr/xmgr file')
        spec.input('tablep_file', valid_type=SinglefileData, required=False, help='xvgr/xmgr file')
        spec.input('rerun_file', valid_type=SinglefileData, required=False, help='Trajectory: xtc trr cpt gro g96 pdb tng')
        spec.input('ei_file', valid_type=SinglefileData, required=False, help='ED sampling input')
        spec.input('multidir_file', valid_type=SinglefileData, required=False, help='Run directory')
        spec.input('awh_file', valid_type=SinglefileData, required=False, help='xvgr/xmgr file')
        spec.input('membed_file', valid_type=SinglefileData, required=False, help='Generic data file')
        spec.input('mp_file', valid_type=SinglefileData, required=False, help='Topology file')
        spec.input('mn_file', valid_type=SinglefileData, required=False, help='Index file')

        # Required outputs.
        spec.output('stdout', valid_type=SinglefileData, help='stdout')
        spec.output('trrfile', valid_type=SinglefileData, help='Output trajectory.')
        spec.output('grofile', valid_type=SinglefileData, help='Output structure file.')
        spec.output('logfile', valid_type=SinglefileData, help='Output log file.')
        spec.output('enfile', valid_type=SinglefileData, help='Output energy file.')

        # Optional outputs.
        spec.output('x_file', required=False, valid_type=SinglefileData, help='Compressed trajectory (tng format or portable xdr format)')
        spec.output('cpo_file', required=False, valid_type=SinglefileData, help='Checkpoint file.')
        spec.output('dhdl_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('field_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('tpi_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('tpid_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('eo_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('px_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('pf_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('ro_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('ra_file', required=False, valid_type=SinglefileData, help='Log file')
        spec.output('rs_file', required=False, valid_type=SinglefileData, help='Log file')
        spec.output('rt_file', required=False, valid_type=SinglefileData, help='Log file')
        spec.output('mtx_file', required=False, valid_type=SinglefileData, help='Hessian Matrix')
        spec.output('if_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')
        spec.output('swap_file', required=False, valid_type=SinglefileData, help='xvgr/xmgr file')

        # Outputs outside of gromacs
        spec.output('logfile_metadata', valid_type=Dict, help='metadata exracted from gromacs logfile')
        #spec.output('test', valid_type=Dict)

        spec.exit_code(300, 'ERROR_MISSING_OUTPUT_FILES', message='Calculation did not produce all expected output files.')

    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files
            needed by the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        codeinfo = CodeInfo()

        # Setup data structures for files.
        input_options = ["tprfile", "cpi_file", "table_file", "tableb_file", "tablep_file", "rerun_file", "ei_file", "multidir_file", "awh_file", "membed_file", "mp_file", "mn_file"]
        output_options = ["c", "e", "g", "o", "x", "cpo", "dhdl", "field", "tpi", "tpid", "eo", "px", "pf", "ro", "ra", "rs", "rt", "mtx", "if", "swap", "logfile_metadata"]
        cmdline_input_files = {}
        input_files = []
        output_files = []

        # Map input files to AiiDA plugin data types.
        for item in input_options:
            if item in self.inputs:
                cmdline_input_files[item] = self.inputs[item].filename
                input_files.append((
                        self.inputs[item].uuid,
                        self.inputs[item].filename,
                        self.inputs[item].filename,
                    ))

        # Add output files to retrieve list.
        output_files.append(self.metadata.options.output_filename)
        for item in output_options:
            if item in self.inputs.parameters:
                output_files.append(self.inputs.parameters[item])

        # Form the commandline.
        codeinfo.cmdline_params = self.inputs.parameters.cmdline_params(cmdline_input_files)
        
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.withmpi = self.inputs.metadata.options.withmpi

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = input_files
        calcinfo.retrieve_list = output_files

        return calcinfo
