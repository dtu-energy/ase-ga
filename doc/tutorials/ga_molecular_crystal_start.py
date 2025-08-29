from ase import Atoms
from ase.data import atomic_numbers
from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import CellBounds, closest_distances_generator

# Number of randomly generated structures
N = 10

# The building blocks
blocks = [('N2', 8)]
# By writing 'N2', the generator will automatically
# get the N2 geometry using ase.build.molecule.

# A guess for the cell volume in Angstrom^3
box_volume = 30.0 * 8

# The cell splitting scheme:
splits = {(2,): 1, (1,): 1}

# The minimal interatomic distances which the
# initial structures must satisfy. We can take these
# a bit larger than usual because these minimal
# distances will only be applied intermolecularly
# (and not intramolecularly):
Z = atomic_numbers['N']
blmin = closest_distances_generator(
    atom_numbers=[Z], ratio_of_covalent_radii=1.3
)

# The bounds for the randomly generated unit cells:
cellbounds = CellBounds(
    bounds={
        'phi': [30, 150],
        'chi': [30, 150],
        'psi': [30, 150],
        'a': [3, 50],
        'b': [3, 50],
        'c': [3, 50],
    }
)

# The familiar 'slab' object, here only providing
# the PBC as there are no atoms or cell vectors
# that need to be applied.
slab = Atoms('', pbc=True)

# create the starting population
sg = StartGenerator(
    slab,
    blocks,
    blmin,
    box_volume=box_volume,
    cellbounds=cellbounds,
    splits=splits,
    number_of_variable_cell_vectors=3,
    test_too_far=False,
)

# Initialize the database
da = PrepareDB(
    db_file_name='gadb.db', simulation_cell=slab, stoichiometry=[Z] * 16
)

# Generate the new structures
# and add them to the database
for i in range(N):
    a = sg.get_new_candidate()
    da.add_unrelaxed_candidate(a)
