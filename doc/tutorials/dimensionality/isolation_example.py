import numpy as np

import ase.build
from ase import Atoms
from ase.geometry.dimensionality import isolate_components
from ase.visualize import view

# build two slabs of different types of MoS2
rep = [4, 4, 1]
a = ase.build.mx2(formula='MoS2', kind='2H', a=3.18, thickness=3.19) * rep
b = ase.build.mx2(formula='MoS2', kind='1T', a=3.18, thickness=3.19) * rep
positions = np.concatenate([a.get_positions(), b.get_positions() + [0, 0, 7]])
numbers = np.concatenate([a.numbers, b.numbers])
cell = a.cell
atoms = Atoms(numbers=numbers, positions=positions, cell=cell, pbc=[1, 1, 1])
atoms.cell[2, 2] = 14.0


# isolate each component in the whole material
result = isolate_components(atoms)
print('counts:', [(k, len(v)) for k, v in sorted(result.items())])

for dim, components in result.items():
    for atoms in components:
        print(dim)
        view(atoms, block=True)
