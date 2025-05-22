from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import write
from ase.optimize import QuasiNewton

# Find the initial and final states for the reaction.

# Set up a (3 x 3) two layer slab of Ru:
a = 2.70
c = 1.59 * a
sqrt3 = 3.0**0.5
bulk = Atoms(
    '2Cu',
    [(0.0, 0.0, 0.0), (1.0 / 3, 1.0 / 3, -0.5 * c)],
    tags=(1, 1),
    pbc=(1, 1, 0),
)
bulk.set_cell([(a, 0, 0), (a / 2, sqrt3 * a / 2, 0), (0, 0, 1)])
slab = bulk.repeat((4, 4, 1))

# Initial state.
# Add the molecule:
x = a / 2.0
y = a * 3.0**0.5 / 6.0
z = 1.8
d = 1.10  # N2 bond length

# Molecular state parallel to the surface:
slab += Atoms('2N', [(x, y, z), (x + sqrt3 * d / 2, y + d / 2, z)])

# Use the EMT calculator for the forces and energies:
slab.calc = EMT()

# We don't want to worry about the Cu degrees of freedom:
mask = [atom.symbol == 'Cu' for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))
relax = QuasiNewton(slab)
relax.run(fmax=0.05)
print('initial state:', slab.get_potential_energy())
write('N2.traj', slab)

# Now the final state.
# Move the second N atom to a neighboring hollow site:
slab[-1].position = (x + a, y, z)
relax.run()
print('final state:  ', slab.get_potential_energy())
write('2N.traj', slab)
