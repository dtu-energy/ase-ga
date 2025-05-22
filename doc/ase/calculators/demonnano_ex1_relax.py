import numpy as np

from ase import Atoms
from ase.calculators.demonnano import DemonNano
from ase.io import write
from ase.optimize import BFGS

d = 0.9775
t = np.pi / 180 * 110.51
mol = Atoms(
    'H2O', positions=[(d, 0, 0), (d * np.cos(t), d * np.sin(t), 0), (0, 0, 0)]
)

input_arguments = {'DFTB': 'SCC', 'CHARGE': '0.0', 'PARAM': 'PTYPE=BIO'}

calc = DemonNano(label='rundir/', input_arguments=input_arguments)
mol.calc = calc

# optimize geometry
dyn = BFGS(mol, trajectory='test.traj')
dyn.run(fmax=0.01)

write('h2o_optimized.xyz', mol)
