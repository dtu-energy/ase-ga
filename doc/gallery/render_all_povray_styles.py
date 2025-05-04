from ase import io
from ase.io.pov import POVRAY
from ase.utils.ptable import ptable

run_povray = True

atoms = ptable(spacing=3)

######
atoms.write('ptable.true_colors.png')
io.write('ptable.xyz', atoms)
#########
styles = list(POVRAY.material_styles_dict)

for style in styles:
    pov_name = f'ptable.{style}.pov'
    ini_name = pov_name.replace('pov', 'ini')

    kwargs = {  # For povray files only
        'textures': len(atoms) * [style],
        'transparent': True,  # Transparent background
        'canvas_width': 1000,  # Width of canvas in pixels
        'camera_type': 'orthographic angle 65',
    }

    generic_projection_settings = {}

    pov_object = io.write(
        pov_name, atoms, **generic_projection_settings, povray_settings=kwargs
    )

    if run_povray:
        pov_object.render()
