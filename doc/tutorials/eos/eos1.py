import numpy as np

from ase import Atoms
from ase.calculators.emt import EMT
from ase.io.trajectory import Trajectory

a = 4.0  # approximate lattice constant
b = a / 2
ag = Atoms(
    'Ag', cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=1, calculator=EMT()
)  # use EMT potential
cell = ag.get_cell()
traj = Trajectory('Ag.traj', 'w')
for x in np.linspace(0.95, 1.05, 5):
    ag.set_cell(cell * x, scale_atoms=True)
    ag.get_potential_energy()
    traj.write(ag)
