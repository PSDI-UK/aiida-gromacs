from aiida.engine import ToContext, WorkChain
from aiida.orm import Code
from aiida.plugins.factories import CalculationFactory, DataFactory

Pdb2gmxCalculation = CalculationFactory('gromacs.pdb2gmx')
EditconfCalculation = CalculationFactory('gromacs.editconf')
SolvateCalculation = CalculationFactory('gromacs.solvate')

Pdb2gmxParameters = DataFactory('gromacs.pdb2gmx')
EditconfParameters = DataFactory('gromacs.editconf')
SolvateParameters = DataFactory('gromacs.solvate')

class SetupWorkChain(WorkChain):
    """WorkChain for setting up a gromacs simulation automatically."""
    
    @classmethod
    def define(cls, spec):
        """Specify workflow recipe."""
        super().define(spec)
        spec.input('code', valid_type=Code)
        spec.input('pdbfile', valid_type=SinglefileData, help='Input structure.')
        spec.input('pdb2gmxparameters', valid_type=Pdb2gmxParameters, help='Command line parameters for gmx pdb2gmx')    
        spec.input('editconfparameters', valid_type=EditconfParameters, help='Command line parameters for gmx editconf')
        spec.input('solvateparameters', valid_type=SolvateParameters, help='Command line parameters for gmx solvate')

        spec.outline(
            cls.pdb2gmx,
            cls.editconf,
            cls.solvate,
        )
    
    
    def pdb2gmx(self):
        """Convert PDB file to forcefield compliant GRO file"""
        inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.pdb2gmxparameters,
            'pdbfile': self.inputs.pdbfile,
            'metadata': {
                'description': 'pdb2gmx job submission with the aiida_gromacs setup workflow',
            },
        }

        future = self.submit(Pdb2gmxCalculation, **inputs)

        return ToContex(pdb2gmx=future)


    def editconf(self):
        """Add simulation box to GRO file"""
        inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.editconfparameters,
            'grofile': self.ctx.pdb2gmx.outputs.outputfile,
            'metadata': {
                'description': 'editconf job submission with the aiida_gromacs setup workflow',
            },
        }

        future = self.submit(EditconfCalculation, **inputs)

        return ToContex(editconf=future)


    def solvate(self):
        """Add solvent to GRO file"""
        inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.solvateparameters,
            'grofile': self.ctx.editconf.outputs.outputfile,
            'topfile': self.ctx.pdb2gmx.outputs.topfile
            'metadata': {
                'description': 'solvate job submission with the aiida_gromacs setup workflow',
            },
        }

        future = self.submit(SolvateCalculation, **inputs)

        return ToContex(solvate=future)
        

