from gpaw import GPAW

from ase import Atoms
from ase.io import Trajectory

atoms = Atoms('N2', positions=[[0, 0, -1], [0, 0, 1]])
atoms.center(vacuum=3.0)

calc = GPAW(mode='lcao', basis='dzp', txt='gpaw.txt')
atoms.calc = calc

traj = Trajectory('binding_curve.traj', 'w')

step = 0.05
nsteps = int(3 / step)

for i in range(nsteps):
    d = 0.5 + i * step
    atoms.positions[1, 2] = atoms.positions[0, 2] + d
    atoms.center(vacuum=3.0)
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    print('distance, energy', d, e)
    print('force', f)
    traj.write(atoms)
