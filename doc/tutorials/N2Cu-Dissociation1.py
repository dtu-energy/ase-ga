from ase import Atoms
from ase.build import add_adsorbate, fcc111
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import write
from ase.optimize import QuasiNewton

# Find the initial and final states for the reaction.

# Set up a (4 x 4) two layer slab of Cu:
slab = fcc111('Cu', size=(4, 4, 2))
slab.set_pbc((1, 1, 0))

# Initial state.
# Add the N2 molecule oriented at 60 degrees:
d = 1.10  # N2 bond length
N2mol = Atoms(
    'N2', positions=[[0.0, 0.0, 0.0], [0.5 * 3**0.5 * d, 0.5 * d, 0.0]]
)
add_adsorbate(slab, N2mol, height=1.0, position='fcc')

# Use the EMT calculator for the forces and energies:
slab.calc = EMT()

# We don't want to worry about the Cu degrees of freedom,
# so fix these atoms:

mask = [atom.symbol == 'Cu' for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))

# Relax the structure
relax = QuasiNewton(slab)
relax.run(fmax=0.05)
print('initial state:', slab.get_potential_energy())
write('N2.traj', slab)

# Now the final state.
# Move the second N atom to a neighboring hollow site:
slab[-1].position[0] = slab[-2].position[0] + 0.25 * slab.cell[0, 0]
slab[-1].position[1] = slab[-2].position[1]
# and relax.
relax.run()
print('final state:  ', slab.get_potential_energy())
write('2N.traj', slab)
