import numpy as np

from ase.build import bulk
from ase.calculators.lj import LennardJones
from ase.constraints import FixSymmetry, UnitCellFilter
from ase.optimize import BFGS
from ase.spacegroup.symmetrize import check_symmetry

# We setup a bcc Al cell - bcc is unstable with LJ potential
# so without constraint this would relax back to an fcc structure
atoms_prim = bulk('Al', 'bcc', a=2 / np.sqrt(3), cubic=True)

# Now we setup a 2x2x2 supercell, and break the symmetry slightly
atoms_init = atoms_prim * [2, 2, 2]
atoms_init.positions[0, 0] += 1.0e-7  # break symmetry by 1e-7

# We use an LJ calculator, and allow the cell and atomic positions to relax
atoms_unsym = atoms_init.copy()
atoms_unsym.calc = LennardJones()
ucf_unsym = UnitCellFilter(atoms_unsym)

dyn = BFGS(ucf_unsym)
print('Initial Energy', atoms_unsym.get_potential_energy())
dyn.run(fmax=0.001)
print('Final Energy', atoms_unsym.get_potential_energy())

# Now we repeat the optimisation with the symmetrization constraint in place
atoms_sym = atoms_init.copy()
atoms_sym.calc = LennardJones()
atoms_sym.set_constraint(FixSymmetry(atoms_sym))
ucf_sym = UnitCellFilter(atoms_sym)

dyn = BFGS(ucf_sym)
print('Initial Energy', atoms_sym.get_potential_energy())
dyn.run(fmax=0.001)
print('Final Energy', atoms_sym.get_potential_energy())

print(
    'position difference',
    np.linalg.norm(atoms_unsym.get_positions() - atoms_sym.get_positions()),
)

# We print out the initial symmetry groups at two different precision levels
print('initial symmetry at precision 1e-6')
check_symmetry(atoms_init, 1.0e-6, verbose=True)
print('initial symmetry at precision 1e-8')
check_symmetry(atoms_init, 1.0e-8, verbose=True)

# Printing the final symmetries shows that
# the "unsym" case relaxes to a lower energy fcc structure
# with a change in spacegroup, while the "sym" case stays as bcc
print('unsym symmetry after relaxation')
d_unsym = check_symmetry(atoms_unsym, verbose=True)
print('sym symmetry after relaxation')
d_sym = check_symmetry(atoms_sym, verbose=True)
