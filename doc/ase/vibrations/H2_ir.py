from gpaw import GPAW, FermiDirac
from gpaw.cluster import Cluster

from ase import optimize
from ase.build import molecule
from ase.vibrations.infrared import Infrared

h = 0.22

atoms = Cluster(molecule('H2'))
atoms.minimal_box(3.5, h=h)

# relax the molecule
calc = GPAW(h=h, occupations=FermiDirac(width=0.1))
atoms.calc = calc
dyn = optimize.FIRE(atoms)
dyn.run(fmax=0.05)
atoms.write('relaxed.traj')

# finite displacement for vibrations
atoms.calc.set(symmetry={'point_group': False})
ir = Infrared(atoms)
ir.run()
