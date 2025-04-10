# creates: WL.png, Ni111slab2x2.png, WL_rot_c.png, WL_rot_a.png, WL_wrap.png, interface-h2o-wrap.png
# flake8: noqa
import runpy

import numpy as np

from ase.build import fcc111
from ase.io import read, write

runpy.run_path('WL.py')

# Use ase.io.read to load atoms object
W = read('WL.traj')
# View the water unit or print the unit cell size.
write('WL.png', W)
# We will need cellW later.
cellW = W.get_cell()
print(cellW)

# We will need as close a lattice match as possible. lets try this slab.
# Using the ase.build module, we make the fcc111 slab.
slab = fcc111('Ni', size=[2, 4, 3], a=3.55, orthogonal=True)
cell = slab.get_cell()
write('Ni111slab2x2.png', slab)
print(cell)

# Rotate the unit cell first to get the close lattice match with the slab.
W.set_cell(
    [[cellW[1, 1], 0, 0], [0, cellW[0, 0], 0], cellW[2]], scale_atoms=False
)
write('WL_rot_c.png', W)

# Now rotate atoms just like the unit cell
W.rotate(90, 'z', center=(0, 0, 0))
write('WL_rot_a.png', W)

# Now we can use wrap
W.wrap()
write('WL_wrap.png', W)

# Match the water lattice to the slab by rescaling
cell1 = np.array([cell[0], cell[1], cellW[2]])
W.set_cell(cell1, scale_atoms=True)
# Set the positions of the water to be 1.5 aangstrom above the slab.
p = slab.get_positions()
W.center(vacuum=p[:, 2].max() + 1.5, axis=2)

# Finally use extend to combine the slab and waterlayer
interface = slab.copy()
interface.extend(W)
interface.center(vacuum=6, axis=2)
write('interface-h2o-wrap.png', interface)
