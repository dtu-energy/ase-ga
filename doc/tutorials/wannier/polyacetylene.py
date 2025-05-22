import numpy as np
from gpaw import GPAW

from ase import Atoms
from ase.dft.kpoints import monkhorst_pack

kpts = monkhorst_pack((13, 1, 1))
calc = GPAW(
    mode='fd',
    h=0.21,
    xc='PBE',
    kpts=kpts,
    nbands=12,
    txt='poly.txt',
    eigensolver='cg',
    convergence={'bands': 9},
    symmetry='off',
)

CC = 1.38
CH = 1.094
a = 2.45
x = a / 2.0
y = np.sqrt(CC**2 - x**2)
atoms = Atoms(
    'C2H2',
    pbc=(True, False, False),
    cell=(a, 8.0, 6.0),
    calculator=calc,
    positions=[[0, 0, 0], [x, y, 0], [x, y + CH, 0], [0, -CH, 0]],
)
atoms.center()
atoms.get_potential_energy()
calc.write('poly.gpw', mode='all')
