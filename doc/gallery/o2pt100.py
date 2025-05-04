# creates: o2pt100.png
import numpy as np

from ase.build import add_adsorbate, fcc100
from ase.io import write

# the metal slab
atoms = fcc100('Pt', size=[4, 10, 3], vacuum=10)
transmittances = [0 for a in atoms]
bonded_atoms = []

upper_layer_idx = [a.index for a in atoms if a.tag == 1]
middle = atoms.positions[upper_layer_idx, :2].max(axis=0) / 2

# the dissociating oxygen... fake some dissociation curve
gas_dist = 1.1
max_height = 8.0
min_height = 1.0
max_dist = 6

# running index for the bonds
index = len(atoms)

for i, x in enumerate(np.linspace(0, 1.5, 6)):
    height = (max_height - min_height) * np.exp(-2 * x) + min_height
    d = np.exp(1.5 * x) / np.exp(1.5**2) * max_dist + gas_dist
    pos = middle + [0, d / 2]
    add_adsorbate(atoms, 'O', height=height, position=pos)
    pos = middle - [0, d / 2]
    add_adsorbate(atoms, 'O', height=height, position=pos)
    transmittances += [x / 2] * 2

    # we want bonds for the first two molecules
    if i < 2:
        bonded_atoms.append([len(atoms) - 1, len(atoms) - 2])

textures = ['ase3' for a in atoms]

# add some semi-transparent bath (only in x/y direction for this example)
cell = atoms.cell

idx = [a.index for a in atoms if a.symbol == 'Pt']

Nbulk = len(idx)
multiples = [0, 1, -1]
for i in multiples:
    for j in multiples:
        if i == j == 0:
            continue
        chunk = atoms[idx]
        chunk.translate(i * cell[0] + j * cell[1])
        atoms += chunk
        transmittances += [0.8] * Nbulk
        textures += ['pale'] * Nbulk

bbox = [-30, 10, 5, 25]

renderer = write(
    'o2pt100.pov',
    atoms,
    rotation='90z,-75x',
    bbox=bbox,
    show_unit_cell=0,
    povray_settings=dict(
        pause=False,
        canvas_width=1024,
        bondatoms=bonded_atoms,
        camera_type='perspective',
        transmittances=transmittances,
        textures=textures,
    ),
)

renderer.render()
