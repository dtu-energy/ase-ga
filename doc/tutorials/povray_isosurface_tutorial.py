from ase.calculators.vasp import VaspChargeDensity
from ase.io import write

spin_cut_off = 0.4
density_cut_off = 0.15

rotation = '24x, 34y, 14z'
# rotation = '0x, 0y, 0z'


vchg = VaspChargeDensity('CHGCAR')
atoms = vchg.atoms[0]

povray_settings = {
    # For povray files only
    'pause': False,  # Pause when done rendering (only if display)
    'transparent': False,  # Transparent background
    'canvas_width': None,  # Width of canvas in pixels
    'canvas_height': 1024,  # Height of canvas in pixels
    'camera_dist': 25.0,  # Distance from camera to front atom
    'camera_type': 'orthographic angle 35',  # 'perspective angle 20'
    'textures': len(atoms) * ['ase3'],
}

# some more options:
# 'image_plane'  : None,  # Distance from front atom to image plane
#                         # (focal depth for perspective)
# 'camera_type'  : 'perspective', # perspective, ultra_wide_angle
# 'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
# 'area_light'   : [(2., 3., 40.) ,# location
#                   'White',       # color
#                   .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
# 'background'   : 'White',        # color
# 'textures'     : tex, # Length of atoms list of texture names
# 'celllinewidth': 0.05, # Radius of the cylinders representing the cell


generic_projection_settings = {
    'rotation': rotation,
    'radii': atoms.positions.shape[0] * [0.3],
    'show_unit_cell': 1,
}

# write returns a renderer object which needs to have the render method called

write(
    'NiO_marching_cubes1.pov',
    atoms,
    **generic_projection_settings,
    povray_settings=povray_settings,
    isosurface_data=dict(density_grid=vchg.chgdiff[0], cut_off=density_cut_off),
).render()

# spin up density, how to specify color and transparency r,g,b,t and a
# material style from the standard ASE set


write(
    'NiO_marching_cubes2.pov',
    atoms,
    **generic_projection_settings,
    povray_settings=povray_settings,
    isosurface_data=dict(
        density_grid=vchg.chgdiff[0],
        cut_off=density_cut_off,
        closed_edges=True,
        color=[0.25, 0.25, 0.80, 0.1],
        material='simple',
    ),
).render()

# spin down density, how to specify a povray material
# that looks like pink jelly
fun_material = """
  material {
    texture {
      pigment { rgbt < 0.8, 0.25, 0.25, 0.5> }
      finish{ diffuse 0.85 ambient 0.99 brilliance 3 specular 0.5 \
roughness 0.001
        reflection { 0.05, 0.98 fresnel on exponent 1.5 }
        conserve_energy
      }
    }
    interior { ior 1.3 }
  }
  photons {
      target
      refraction on
      reflection on
      collect on
  }"""

write(
    'NiO_marching_cubes3.pov',
    atoms,
    **generic_projection_settings,
    povray_settings=povray_settings,
    isosurface_data=dict(
        density_grid=vchg.chgdiff[0],
        cut_off=-spin_cut_off,
        gradient_ascending=True,
        material=fun_material,
    ),
).render()
