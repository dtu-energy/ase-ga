import sys

from ase.calculators.emt import EMT
from ase.ga.relax_attaches import VariansBreak
from ase.io import read, write
from ase.optimize import BFGS

fname = sys.argv[1]

print(f'Now relaxing {fname}')
a = read(fname)

a.calc = EMT()
dyn = BFGS(a, trajectory=None, logfile=None)
vb = VariansBreak(a, dyn)
dyn.attach(vb.write)
dyn.run(fmax=0.05)

a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()

write(fname[:-5] + '_done.traj', a)

print(f'Done relaxing {fname}')
