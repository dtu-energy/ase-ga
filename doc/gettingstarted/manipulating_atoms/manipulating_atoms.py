# creates: a1.png a2.png a3.png
# creates: WL.png, Ni111slab2x2.png, WL_rot_c.png, WL_rot_a.png, WL_wrap.png
# creates: interface-h2o-wrap.png

import runpy
from math import sqrt

import numpy as np

from ase import Atoms
from ase.build import fcc111
from ase.io import read, write

a = 3.55
atoms = Atoms(
    'Ni4',
    cell=[sqrt(2) * a, sqrt(2) * a, 1.0, 90, 90, 120],
    pbc=(1, 1, 0),
    scaled_positions=[(0, 0, 0), (0.5, 0, 0), (0, 0.5, 0), (0.5, 0.5, 0)],
)
atoms.center(vacuum=5.0, axis=2)
write('a1.png', atoms, rotation='-73x')
a3 = atoms.repeat((3, 3, 2))
a3.cell = atoms.cell
write('a2.png', a3, rotation='-73x')
h = 1.9
relative = (1 / 6, 1 / 6, 0.5)
absolute = np.dot(relative, atoms.cell) + (0, 0, h)
atoms.append('Ag')
atoms.positions[-1] = absolute
a3 = atoms.repeat((3, 3, 2))
a3.cell = atoms.cell
write('a3.png', a3)

runpy.run_path('WL.py')
W = read('WL.traj')
write('WL.png', W)
slab = fcc111('Ni', size=[2, 4, 3], a=3.55, orthogonal=True)
write('Ni111slab2x2.png', slab)
W.cell = [W.cell[1, 1], W.cell[0, 0], 0.0]
write('WL_rot_c.png', W)
W.rotate(90, 'z', center=(0, 0, 0))
write('WL_rot_a.png', W)
W.wrap()
write('WL_wrap.png', W)
W.set_cell(slab.cell, scale_atoms=True)
zmin = W.positions[:, 2].min()
zmax = slab.positions[:, 2].max()
W.positions += (0, 0, zmax - zmin + 1.5)
interface = slab + W
interface.center(vacuum=6, axis=2)
write('interface-h2o-wrap.png', interface)
