import numpy as np
from gpaw import GPAW, PW

from ase.build import bulk
from ase.calculators.dftd3 import DFTD3
from ase.constraints import UnitCellFilter
from ase.optimize import LBFGS

rng = np.random.RandomState(0)

diamond = bulk('C')
diamond.rattle(stdev=0.1, rng=rng)
diamond.cell += rng.normal(scale=0.1, size=(3, 3))
dft = GPAW(xc='PBE', kpts=(8, 8, 8), mode=PW(400))
d3 = DFTD3(dft=dft)
diamond.calc = d3

ucf = UnitCellFilter(diamond)

opt = LBFGS(ucf, logfile='diamond_opt.log', trajectory='diamond_opt.traj')
opt.run(fmax=0.05)
