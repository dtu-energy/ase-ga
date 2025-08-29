from ase import Atoms
from ase.data import atomic_numbers
from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import CellBounds, closest_distances_generator

# Number of random initial structures to generate
N = 20

# Target cell volume for the initial structures, in angstrom^3
volume = 240.0

# Specify the 'building blocks' from which the initial structures
# will be constructed. Here we take single Ag atoms as building
# blocks, 24 in total.
blocks = [('Ag', 24)]
# We may also write:
blocks = ['Ag'] * 24

# Generate a dictionary with the closest allowed interatomic distances
Z = atomic_numbers['Ag']
blmin = closest_distances_generator(
    atom_numbers=[Z], ratio_of_covalent_radii=0.5
)

# Specify reasonable bounds on the minimal and maximal
# cell vector lengths (in angstrom) and angles (in degrees)
cellbounds = CellBounds(
    bounds={
        'phi': [35, 145],
        'chi': [35, 145],
        'psi': [35, 145],
        'a': [3, 50],
        'b': [3, 50],
        'c': [3, 50],
    }
)

# Choose an (optional) 'cell splitting' scheme which basically
# controls the level of translational symmetry (within the unit
# cell) of the randomly generated structures. Here a 1:1 ratio
# of splitting factors 2 and 1 is used:
splits = {(2,): 1, (1,): 1}
# There will hence be a 50% probability that a candidate
# is constructed by repeating an randomly generated Ag12
# structure along a randomly chosen axis. In the other 50%
# of cases, no cell cell splitting will be applied.

# The 'slab' object in the GA serves as a template
# in the creation of new structures, which inherit
# the slab's atomic positions (if any), cell vectors
# (if specified), and periodic boundary conditions.
# Here only the last property is relevant:
slab = Atoms('', pbc=True)

# Initialize the random structure generator
sg = StartGenerator(
    slab,
    blocks,
    blmin,
    box_volume=volume,
    number_of_variable_cell_vectors=3,
    cellbounds=cellbounds,
    splits=splits,
)

# Create the database
da = PrepareDB(db_file_name='gadb.db', stoichiometry=[Z] * 24)

# Generate N random structures
# and add them to the database
for i in range(N):
    a = sg.get_new_candidate()
    da.add_unrelaxed_candidate(a)
