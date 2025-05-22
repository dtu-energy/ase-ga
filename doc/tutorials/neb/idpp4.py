import numpy as np

from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.lattice.cubic import FaceCenteredCubic
from ase.mep import NEB
from ase.optimize.fire import FIRE as QuasiNewton

# Set the number of images you want.
nimages = 5

# Some algebra to determine surface normal and the plane of the surface.
d3 = [2.0, 1.0, 1.0]
a1 = np.array([0.0, 1.0, 1.0])
d1 = np.cross(a1, d3)
a2 = np.array([0.0, -1.0, 1.0])
d2 = np.cross(a2, d3)

# Create the slab.
slab = FaceCenteredCubic(
    directions=[d1, d2, d3], size=(2, 1, 2), symbol=('Pt'), latticeconstant=3.9
)

# Add some vacuum to the slab.
uc = slab.get_cell()
uc[2] += [0.0, 0.0, 10.0]  # There are ten layers of vacuum.
uc = slab.set_cell(uc, scale_atoms=False)

# Some positions needed to place the atom in the correct place.
x1 = 1.379
x2 = 4.137
x3 = 2.759
y1 = 0.0
y2 = 2.238
z1 = 7.165
z2 = 6.439


# Add the adatom to the list of atoms and set constraints of surface atoms.
slab += Atoms('N', [((x2 + x1) / 2, y1, z1 + 1.5)])
mask = [atom.symbol == 'Pt' for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))

# Optimise the initial state: atom below step.
initial = slab.copy()
initial.calc = EMT()
relax = QuasiNewton(initial)
relax.run(fmax=0.05)

# Optimise the final state: atom above step.
slab[-1].position = (x3, y2 + 1.0, z2 + 3.5)
final = slab.copy()
final.calc = EMT()
relax = QuasiNewton(final)
relax.run(fmax=0.05)

# Create a list of images for interpolation.
images = [initial]
for i in range(nimages):
    images.append(initial.copy())

for image in images:
    image.calc = EMT()

images.append(final)

# Carry out linear interpolation.
neb = NEB(images)
neb.interpolate()

# Run NEB calculation.
qn = QuasiNewton(
    neb, trajectory='N_diffusion_lin.traj', logfile='N_diffusion_lin.log'
)
qn.run(fmax=0.05)
