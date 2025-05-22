from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import read
from ase.mep import NEB
from ase.optimize import BFGS
from ase.parallel import world

initial = read('initial.traj')
final = read('final.traj')

constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])

images = [initial]
j = world.rank * 3 // world.size  # my image number
for i in range(3):
    image = initial.copy()
    if i == j:
        image.calc = EMT()
    image.set_constraint(constraint)
    images.append(image)
images.append(final)

neb = NEB(images, parallel=True)
neb.interpolate()
qn = BFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)
