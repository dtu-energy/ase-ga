# creates: precon.png

import matplotlib.pyplot as plt

from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.loggingcalc import LoggingCalculator
from ase.optimize.precon import Exp, PreconLBFGS

a0 = bulk('Cu', cubic=True)
a0 *= [3, 3, 3]
del a0[0]
a0.rattle(0.1)

nsteps = []
energies = []
log_calc = LoggingCalculator(EMT())

for precon, label in [(None, 'None'), (Exp(A=3, mu=1.0), 'Exp(A=3)')]:
    log_calc.label = label
    atoms = a0.copy()
    atoms.calc = log_calc
    opt = PreconLBFGS(atoms, precon=precon, use_armijo=True)
    opt.run(fmax=1e-3)

log_calc.plot(markers=['r-', 'b-'], energy=False, lw=2)
plt.savefig('precon.png')
