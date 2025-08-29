from ga_bulk_relax import relax

from ase.ga import get_raw_score
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.ofp_comparator import OFPComparator
from ase.ga.population import Population
from ase.ga.soft_mutation import SoftMutation
from ase.ga.standardmutations import StrainMutation
from ase.ga.utilities import CellBounds, closest_distances_generator
from ase.io import write

# Connect to the database and retrieve some information
da = DataConnection('gadb.db')
slab = da.get_slab()
atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_top = len(atom_numbers_to_optimize)

# Use Oganov's fingerprint functions to decide whether
# two structures are identical or not
comp = OFPComparator(
    n_top=n_top,
    dE=1.0,
    cos_dist_max=1e-3,
    rcut=10.0,
    binwidth=0.05,
    pbc=[True, True, True],
    sigma=0.05,
    nsigma=4,
    recalculate=False,
)

# Define the cell and interatomic distance bounds
# that the candidates must obey
blmin = closest_distances_generator(atom_numbers_to_optimize, 0.5)

cellbounds = CellBounds(
    bounds={
        'phi': [20, 160],
        'chi': [20, 160],
        'psi': [20, 160],
        'a': [2, 60],
        'b': [2, 60],
        'c': [2, 60],
    }
)

# Define a pairing operator with 100% (0%) chance that the first
# (second) parent will be randomly translated, and with each parent
# contributing to at least 15% of the child's scaled coordinates
pairing = CutAndSplicePairing(
    slab,
    n_top,
    blmin,
    p1=1.0,
    p2=0.0,
    minfrac=0.15,
    number_of_variable_cell_vectors=3,
    cellbounds=cellbounds,
    use_tags=False,
)

# Define a strain mutation with a typical standard deviation of 0.7
# for the strain matrix elements (drawn from a normal distribution)
strainmut = StrainMutation(
    blmin,
    stddev=0.7,
    cellbounds=cellbounds,
    number_of_variable_cell_vectors=3,
    use_tags=False,
)

# Define a soft mutation; we need to provide a dictionary with
# (typically rather short) minimal interatomic distances which
# is used to determine when to stop displacing the atoms along
# the chosen mode. The minimal and maximal single-atom displacement
# distances (in Angstrom) for a valid mutation are provided via
# the 'bounds' keyword argument.
blmin_soft = closest_distances_generator(atom_numbers_to_optimize, 0.1)
softmut = SoftMutation(blmin_soft, bounds=[2.0, 5.0], use_tags=False)
# By default, the operator will update a "used_modes.json" file
# after every mutation, listing which modes have been used so far
# for each structure in the database. The mode indices start at 3
# as the three lowest frequency modes are translational modes.

# Set up the relative probabilities for the different operators
operators = OperationSelector([4.0, 3.0, 3.0], [pairing, softmut, strainmut])

# Relax the initial candidates
while da.get_number_of_unrelaxed_candidates() > 0:
    a = da.get_an_unrelaxed_candidate()

    relax(a, cellbounds=cellbounds)
    da.add_relaxed_step(a)

    cell = a.get_cell()
    if not cellbounds.is_within_bounds(cell):
        da.kill_candidate(a.info['confid'])

# Initialize the population
population_size = 20
population = Population(
    data_connection=da,
    population_size=population_size,
    comparator=comp,
    logfile='log.txt',
    use_extinct=True,
)

# Update the scaling volume used in some operators
# based on a number of the best candidates
current_pop = population.get_current_population()
strainmut.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)
pairing.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)

# Test n_to_test new candidates; in this example we need
# only few GA iterations as the global minimum (FCC Ag)
# is very easily found (typically already after relaxation
# of the initial random structures).
n_to_test = 50

for step in range(n_to_test):
    print(f'Now starting configuration number {step}')

    # Create a new candidate
    a3 = None
    while a3 is None:
        a1, a2 = population.get_two_candidates()
        a3, desc = operators.get_new_individual([a1, a2])

    # Save the unrelaxed candidate
    da.add_unrelaxed_candidate(a3, description=desc)

    # Relax the new candidate and save it
    relax(a3, cellbounds=cellbounds)
    da.add_relaxed_step(a3)

    # If the relaxation has changed the cell parameters
    # beyond the bounds we disregard it in the population
    cell = a3.get_cell()
    if not cellbounds.is_within_bounds(cell):
        da.kill_candidate(a3.info['confid'])

    # Update the population
    population.update()

    if step % 10 == 0:
        # Update the scaling volumes of the strain mutation
        # and the pairing operator based on the current
        # best structures contained in the population
        current_pop = population.get_current_population()
        strainmut.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)
        pairing.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)
        write('current_population.traj', current_pop)

print('GA finished after step %d' % step)
hiscore = get_raw_score(current_pop[0])
print('Highest raw score = %8.4f eV' % hiscore)

all_candidates = da.get_all_relaxed_candidates()
write('all_candidates.traj', all_candidates)

current_pop = population.get_current_population()
write('current_population.traj', current_pop)
