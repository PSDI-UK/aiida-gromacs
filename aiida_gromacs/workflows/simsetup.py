from aiida.engine import ToContext, WorkChain
from aiida.orm import Code, SinglefileData
from aiida.plugins.factories import CalculationFactory, DataFactory

Pdb2gmxCalculation = CalculationFactory('gromacs.pdb2gmx')
EditconfCalculation = CalculationFactory('gromacs.editconf')
SolvateCalculation = CalculationFactory('gromacs.solvate')
GromppCalculation = CalculationFactory('gromacs.grompp')
GenionCalculation = CalculationFactory('gromacs.genion')

Pdb2gmxParameters = DataFactory('gromacs.pdb2gmx')
EditconfParameters = DataFactory('gromacs.editconf')
SolvateParameters = DataFactory('gromacs.solvate')
GromppParameters = DataFactory('gromacs.grompp')
GenionParameters = DataFactory('gromacs.genion')

class SetupWorkChain(WorkChain):
    """WorkChain for setting up a gromacs simulation automatically."""
    
    @classmethod
    def define(cls, spec):
        """Specify workflow recipe."""
        super().define(spec)
        spec.input('code', valid_type=Code)
        spec.input('pdbfile', valid_type=SinglefileData, help='Input structure.')
        spec.input('ionsmdp', valid_type=SinglefileData, help='Input structure.')
        spec.input('pdb2gmxparameters', valid_type=Pdb2gmxParameters, help='Command line parameters for gmx pdb2gmx')    
        spec.input('editconfparameters', valid_type=EditconfParameters, help='Command line parameters for gmx editconf')
        spec.input('solvateparameters', valid_type=SolvateParameters, help='Command line parameters for gmx solvate')
        spec.input('gromppionsparameters', valid_type=GromppParameters, help='Command line parameters for gmx grompp')
        spec.input('genionparameters', valid_type=GenionParameters, help='Command line parameters for gmx genion')

        spec.outline(
            cls.pdb2gmx,
            cls.editconf,
            cls.solvate,
            cls.gromppions,
            cls.genion,
            cls.result,
        )
        
        spec.output('result')
    
    
    def pdb2gmx(self):
        """Convert PDB file to forcefield compliant GRO file"""
        inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.pdb2gmxparameters,
            'pdbfile': self.inputs.pdbfile,
            'metadata': {
                'description': 'convert pdb file to gromacs gro format.',
            },
        }

        future = self.submit(Pdb2gmxCalculation, **inputs)

        return ToContext(pdb2gmx=future)


    def editconf(self):
        """Add simulation box to GRO file."""
        inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.editconfparameters,
            'grofile': self.ctx.pdb2gmx.outputs.outputfile,
            'metadata': {
                'description': 'add simulation box parameters to gro file.',
            },
        }

        future = self.submit(EditconfCalculation, **inputs)

        return ToContext(editconf=future)


    def solvate(self):
        """Add solvent to GRO file."""
        inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.solvateparameters,
            'grofile': self.ctx.editconf.outputs.outputfile,
            'topfile': self.ctx.pdb2gmx.outputs.topfile,
            'metadata': {
                'description': 'add solvent to simulation box.',
            },
        }

        future = self.submit(SolvateCalculation, **inputs)

        return ToContext(solvate=future)


    def gromppions(self):
        """Create a tpr for adding ions."""
        inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.gromppionsparameters,
            'mdpfile': self.inputs.ionsmdp,
            'grofile': self.ctx.solvate.outputs.outputfile,
            'topfile': self.ctx.solvate.outputs.topfile,
            'metadata': {
                'description': 'prepare the tpr for adding ions.',
            },
        }

        future = self.submit(GromppCalculation, **inputs)

        return ToContext(grommp=future)


    def genion(self):
        """Add ions to system to balance charge."""
        
        #TODO sort this out.
        from aiida_gromacs import helpers
        computer = helpers.get_computer()
        gromacs_code = helpers.get_code(entry_point='bash',
                                        computer=computer)
        inputs = {
            'code': gromacs_code,
            'parameters': self.inputs.genionparameters,
            'tprfile': self.ctx.grommp.outputs.outputfile,
            'topfile': self.ctx.solvate.outputs.topfile,
            'metadata': {
                'description': 'add ions to simulation box.',
            },
        }

        future = self.submit(GenionCalculation, **inputs)

        return ToContext(genion=future)


    def result(self):
        """Results"""
        self.out('result', self.ctx.genion.outputs.outputfile)
