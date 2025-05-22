from ase.calculators.emt import EMT
from ase.collections import dcdft
from ase.io import Trajectory

for symbol in ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']:
    traj = Trajectory(f'{symbol}.traj', 'w')
    for s in range(94, 108, 2):
        atoms = dcdft[symbol]
        atoms.set_cell(atoms.cell * (s / 100) ** (1 / 3), scale_atoms=True)
        atoms.calc = EMT()
        atoms.get_potential_energy()
        traj.write(atoms)
