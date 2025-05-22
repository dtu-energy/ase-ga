# creates: atomic_radii.png

import matplotlib.pyplot as plt
import numpy as np

from ase.data import chemical_symbols, covalent_radii
from ase.data.vdw import vdw_radii as vdw1
from ase.data.vdw_alvarez import vdw_radii as vdw2

plt.grid(ls=':')
c1 = covalent_radii.copy()
c1[c1 < 0.2001] = np.nan  # Remove 'false' values which are all 0.2
plt.plot(vdw2, marker='.', label='vdw_radii [ase.data.vdw_alvarez]')
plt.plot(vdw1, marker='.', label='vdw_radii [ase.data.vdw]')
plt.plot(c1, marker='.', label='covalent_radii [ase.data]')
nobles = [2, 10, 18, 36, 54, 86]
plt.xticks(nobles, [chemical_symbols[Z] for Z in nobles])
plt.xlabel('Z')
plt.ylabel('radius [Å]')
plt.legend(loc='best')
plt.savefig('atomic_radii.png')
