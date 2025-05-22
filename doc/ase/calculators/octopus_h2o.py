from ase.calculators.octopus import Octopus
from ase.collections import g2
from ase.optimize import QuasiNewton

# Water molecule with somewhat randomized initial positions:
atoms = g2['H2O']
atoms.rattle(stdev=0.1, seed=42)

calc = Octopus(directory='oct-h2o', Spacing=0.25)
atoms.calc = calc

opt = QuasiNewton(atoms, logfile='opt.log', trajectory='opt.traj')
opt.run(fmax=0.05)
