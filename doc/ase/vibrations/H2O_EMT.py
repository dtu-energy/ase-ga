from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.vibrations import Vibrations

h2o = molecule('H2O')
h2o.calc = EMT()
BFGS(h2o).run(fmax=0.01)

vib = Vibrations(h2o)
vib.run()
vib.summary(log='H2O_EMT_summary.txt')
vib.write_mode(-1)
