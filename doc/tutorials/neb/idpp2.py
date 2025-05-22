from ase.build import molecule
from ase.calculators.emt import EMT
from ase.mep import NEB
from ase.optimize.fire import FIRE as QuasiNewton

# Optimise molecule.
initial = molecule('C2H6')
initial.calc = EMT()
relax = QuasiNewton(initial)
relax.run(fmax=0.05)

# Create final state.
final = initial.copy()
final.positions[2:5] = initial.positions[[3, 4, 2]]

# Generate blank images.
images = [initial]

for i in range(9):
    images.append(initial.copy())

for image in images:
    image.calc = EMT()

images.append(final)

# Run IDPP interpolation.
neb = NEB(images)
neb.interpolate('idpp')

# Run NEB calculation.
qn = QuasiNewton(neb, trajectory='ethane_idpp.traj', logfile='ethane_idpp.log')
qn.run(fmax=0.05)
