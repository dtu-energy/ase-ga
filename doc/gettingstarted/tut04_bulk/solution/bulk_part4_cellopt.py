from gpaw import GPAW, PW

from ase.constraints import ExpCellFilter
from ase.io import write
from ase.optimize import BFGS
from ase.spacegroup import crystal

a = 4.6
c = 2.95

# Rutile TiO2:
atoms = crystal(
    ['Ti', 'O'],
    basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
    spacegroup=136,
    cellpar=[a, a, c, 90, 90, 90],
)
write('rutile.traj', atoms)

calc = GPAW(mode=PW(800), kpts=[2, 2, 3], txt='gpaw.rutile.txt')
atoms.calc = calc

opt = BFGS(ExpCellFilter(atoms), trajectory='opt.rutile.traj')
opt.run(fmax=0.05)

calc.write('groundstate.rutile.gpw')

print('Final lattice:')
print(atoms.cell.get_bravais_lattice())
