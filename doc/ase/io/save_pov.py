# creates: NaCl_C6H6.png

import numpy as np

from ase import Atoms
from ase.build import molecule
from ase.io import write

a = 5.64  # Lattice constant for NaCl
cell = [a / np.sqrt(2), a / np.sqrt(2), a]
atoms = Atoms(
    symbols='Na2Cl2',
    pbc=True,
    cell=cell,
    scaled_positions=[
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.5),
        (0.5, 0.5, 0.0),
        (0.0, 0.0, 0.5),
    ],
) * (3, 4, 2) + molecule('C6H6')

# Move molecule to 3.5Ang from surface, and translate one unit cell in xy
atoms.positions[-12:, 2] += atoms.positions[:-12, 2].max() + 3.5
atoms.positions[-12:, :2] += cell[:2]

# Mark a single unit cell
atoms.cell = cell

# View used to start ag, and find desired viewing angle
# view(atoms)
rot = '35x,63y,36z'  # found using ag: 'view -> rotate'

# Common kwargs for eps, png, pov
generic_projection_settings = {
    'rotation': rot,  # text string with rotation (default='' )
    'radii': 0.85,  # float, or a list with one float per atom
    'colors': None,  # List: one (r, g, b) tuple per atom
    'show_unit_cell': 2,  # 0, 1, or 2 to not show, show, and show all of cell
}

# Extra kwargs only available for povray (All units in angstrom)
povray_settings = {
    'display': False,  # Display while rendering
    'pause': True,  # Pause when done rendering (only if display)
    'transparent': False,  # Transparent background
    'canvas_width': None,  # Width of canvas in pixels
    'canvas_height': None,  # Height of canvas in pixels
    'camera_dist': 50.0,  # Distance from camera to front atom
    'image_plane': None,  # Distance from front atom to image plane
    'camera_type': 'perspective',  # perspective, ultra_wide_angle
    'point_lights': [],  # [[loc1, color1], [loc2, color2],...]
    'area_light': [
        (2.0, 3.0, 40.0),  # location
        'White',  # color
        0.7,
        0.7,
        3,
        3,
    ],  # width, height, Nlamps_x, Nlamps_y
    'background': 'White',  # color
    'textures': None,  # Length of atoms list of texture names
    'celllinewidth': 0.1,  # Radius of the cylinders representing the cell
}

# Write the .pov (and .ini) file.
# comment out render not call the povray executable
renderer = write(
    'NaCl_C6H6.pov',
    atoms,
    **generic_projection_settings,
    povray_settings=povray_settings,
)

renderer.render()
