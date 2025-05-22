from ase import Atoms
from ase.calculators.exciting.exciting import ExcitingGroundStateTemplate

# test structure, not real
nitrogen_trioxide_atoms = Atoms(
    'NO3',
    cell=[[2, 2, 0], [0, 4, 0], [0, 0, 6]],
    positions=[(0, 0, 0), (1, 3, 0), (0, 0, 1), (0.5, 0.5, 0.5)],
    pbc=True,
)

gs_template_obj = ExcitingGroundStateTemplate()

# Write an exciting input.xml file for the NO3 system.
gs_template_obj.write_input(
    directory='./',
    atoms=nitrogen_trioxide_atoms,
    parameters={
        'title': None,
        'species_path': './',
        'ground_state_input': {
            'rgkmax': 8.0,
            'do': 'fromscratch',
            'ngridk': [6, 6, 6],
            'xctype': 'GGA_PBE_SOL',
            'vkloff': [0, 0, 0],
        },
    },
)
