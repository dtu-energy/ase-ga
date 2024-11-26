import numpy as np
from ga_molecular_crystal_relax import relax

from ase.ga import get_raw_score
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.ofp_comparator import OFPComparator
from ase.ga.population import Population
from ase.ga.soft_mutation import SoftMutation
from ase.ga.standardmutations import (
    RattleMutation,
    RattleRotationalMutation,
    RotationalMutation,
    StrainMutation,
)
from ase.ga.utilities import CellBounds, closest_distances_generator
from ase.io import write

da = DataConnection('gadb.db')

# Various items needed for initializing the genetic operators
slab = da.get_slab()
atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_top = len(atom_numbers_to_optimize)
blmin = closest_distances_generator(atom_numbers_to_optimize, 1.0)
cellbounds = CellBounds(bounds={'phi': [30, 150], 'chi': [30, 150],
                                'psi': [30, 150]})

# Note the "use_tags" keyword argument being used
# to signal that we want to preserve molecular identity
# via the tags
pairing = CutAndSplicePairing(slab, n_top, blmin, p1=1., p2=0.,
                              minfrac=0.15, cellbounds=cellbounds,
                              number_of_variable_cell_vectors=3,
                              use_tags=True)

rattlemut = RattleMutation(blmin, n_top, rattle_prop=0.3, rattle_strength=0.5,
                           use_tags=True)

strainmut = StrainMutation(blmin, stddev=0.7, cellbounds=cellbounds,
                           use_tags=True)

rotmut = RotationalMutation(blmin, fraction=0.3, min_angle=0.5 * np.pi)

rattlerotmut = RattleRotationalMutation(rattlemut, rotmut)

blmin_soft = closest_distances_generator(atom_numbers_to_optimize, 0.8)
softmut = SoftMutation(blmin_soft, bounds=[2., 5.], use_tags=True)

operators = OperationSelector([5, 1, 1, 1, 1, 1], [
    pairing, rattlemut, strainmut, rotmut, rattlerotmut, softmut])

# Relaxing the initial candidates
while da.get_number_of_unrelaxed_candidates() > 0:
    a = da.get_an_unrelaxed_candidate()
    relax(a)
    da.add_relaxed_step(a)

# The structure comparator for the population
comp = OFPComparator(n_top=n_top, dE=1.0, cos_dist_max=5e-3, rcut=10.,
                     binwidth=0.05, pbc=[True, True, True], sigma=0.05,
                     nsigma=4, recalculate=False)

# The population
population = Population(data_connection=da,
                        population_size=10,
                        comparator=comp,
                        logfile='log.txt')

current_pop = population.get_current_population()
strainmut.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)
pairing.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)

# Test a few new candidates
n_to_test = 10

for step in range(n_to_test):
    print(f'Now starting configuration number {step}', flush=True)

    # Generate a new candidate
    a3 = None
    while a3 is None:
        a1, a2 = population.get_two_candidates()
        a3, desc = operators.get_new_individual([a1, a2])

    # Relax it and add to database
    da.add_unrelaxed_candidate(a3, description=desc)
    relax(a3)
    da.add_relaxed_step(a3)

    # Update the population
    population.update()
    current_pop = population.get_current_population()
    write('current_population.traj', current_pop)

    # Update the strain mutation and pairing operators
    if step % 10 == 0:
        strainmut.update_scaling_volume(current_pop, w_adapt=0.5,
                                        n_adapt=4)
        pairing.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)

    # Print out information for easier follow-up/analysis/plotting:
    print('Step %d %s %.3f %.3f %.3f' % (
        step, desc, get_raw_score(a1), get_raw_score(a2), get_raw_score(a3)))

    print('Step %d highest raw score in pop: %.3f' %
          (step, get_raw_score(current_pop[0])))

print('GA finished after step %d' % step)
hiscore = get_raw_score(current_pop[0])
print('Highest raw score = %8.4f eV' % hiscore)
write('all_candidates.traj', da.get_all_relaxed_candidates())
write('current_population.traj', current_pop)
