from aiida.engine import WorkChain
from aiida.plugins.factories import CalculationFactory

Pdb2gmxCalculation = CalculationFactory('gromacs.pdb2gmx')
EditconfCalculation = CalculationFactory('gromacs.editconf')
SolvateCalculation = CalculationFactory('gromacs.solvate')

class SetupWorkChain(WorkChain):
    """WorkChain for setting up a gromacs simulation automatically."""
    
    @classmethod
    def define(cls, spec):
    """Specify workflow recipe."""
    super().define(spec)
    
    # TODO input spec
    
    
    spec.outline(
        cls.pdb2gmx,
        cls.editconf,
        cls.solvate,
        )
    
    
    def pdb2gmx(self):
        pass
    

    def editconf(self):
        pass


    def solvate(self):
        pass
