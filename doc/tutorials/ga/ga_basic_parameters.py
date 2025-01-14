from random import random

import numpy as np

from ase.ga import get_parametrization
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.pbs_queue_run import PBSQueueRun
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.standardmutations import (
    MirrorMutation,
    PermutationMutation,
    RattleMutation,
)
from ase.ga.utilities import (
    closest_distances_generator,
    get_all_atom_types,
    get_angles_distribution,
    get_atoms_connections,
    get_atoms_distribution,
    get_neighborlist,
    get_rings,
)
from ase.io import write


def jtg(job_name, traj_file):
    s = '#!/bin/sh\n'
    s += '#PBS -l nodes=1:ppn=16\n'
    s += '#PBS -l walltime=100:00:00\n'
    s += f'#PBS -N {job_name}\n'
    s += '#PBS -q q16\n'
    s += 'cd $PBS_O_WORKDIR\n'
    s += 'NPROCS==`wc -l < $PBS_NODEFILE`\n'
    s += 'mpirun --mca mpi_warn_on_fork 0 -np $NPROCS '
    s += f'gpaw-python calc_gpaw.py {traj_file}\n'
    return s


def combine_parameters(conf):
    # Get and combine selected parameters
    parameters = []
    gets = [get_atoms_connections(conf) + get_rings(conf) +
            get_angles_distribution(conf) + get_atoms_distribution(conf)]
    for get in gets:
        parameters += get
    return parameters


def should_we_skip(conf, comparison_energy, weights):
    parameters = combine_parameters(conf)
    # Return if weights not defined (too few completed
    # calculated structures to make a good fit)
    if weights is None:
        return False
    regression_energy = sum(p * q for p, q in zip(weights, parameters))
    # Skip with 90% likelihood if energy appears to go up 5 eV or more
    if (regression_energy - comparison_energy) > 5 and random() < 0.9:
        return True
    else:
        return False


population_size = 20
mutation_probability = 0.3

# Initialize the different components of the GA
da = DataConnection('gadb.db')
tmp_folder = 'work_folder/'
# The PBS queing interface is created
pbs_run = PBSQueueRun(da,
                      tmp_folder=tmp_folder,
                      job_prefix='Ag2Au2_opt',
                      n_simul=5,
                      job_template_generator=jtg,
                      find_neighbors=get_neighborlist,
                      perform_parametrization=combine_parameters)

atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
slab = da.get_slab()
all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
blmin = closest_distances_generator(all_atom_types,
                                    ratio_of_covalent_radii=0.7)

comp = InteratomicDistanceComparator(n_top=n_to_optimize,
                                     pair_cor_cum_diff=0.015,
                                     pair_cor_max=0.7,
                                     dE=0.02,
                                     mic=False)
pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
mutations = OperationSelector([1., 1., 1.],
                              [MirrorMutation(blmin, n_to_optimize),
                               RattleMutation(blmin, n_to_optimize),
                               PermutationMutation(n_to_optimize)])

# Relax all unrelaxed structures (e.g. the starting population)
while (da.get_number_of_unrelaxed_candidates() > 0 and
       not pbs_run.enough_jobs_running()):
    a = da.get_an_unrelaxed_candidate()
    pbs_run.relax(a)


# create the population
population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

# create the regression expression for estimating the energy
all_trajs = da.get_all_relaxed_candidates()
sampled_points = []
sampled_energies = []
for conf in all_trajs:
    no_of_conn = list(get_parametrization(conf))
    if no_of_conn not in sampled_points:
        sampled_points.append(no_of_conn)
        sampled_energies.append(conf.get_potential_energy())

sampled_points = np.array(sampled_points)
sampled_energies = np.array(sampled_energies)

if len(sampled_points) > 0 and len(sampled_energies) >= len(sampled_points[0]):
    weights = np.linalg.lstsq(sampled_points, sampled_energies, rcond=-1)[0]
else:
    weights = None

# Submit new candidates until enough are running
while (not pbs_run.enough_jobs_running() and
       len(population.get_current_population()) > 2):
    a1, a2 = population.get_two_candidates()

    # Selecting the "worst" parent energy
    # which the child should be compared to
    ce_a1 = da.get_atoms(a1.info['relax_id']).get_potential_energy()
    ce_a2 = da.get_atoms(a2.info['relax_id']).get_potential_energy()
    comparison_energy = min(ce_a1, ce_a2)

    a3, desc = pairing.get_new_individual([a1, a2])
    if a3 is None:
        continue
    if should_we_skip(a3, comparison_energy, weights):
        continue
    da.add_unrelaxed_candidate(a3, description=desc)

    if random() < mutation_probability:
        a3_mut, desc_mut = mutations.get_new_individual([a3])
        if (a3_mut is not None and
                not should_we_skip(a3_mut, comparison_energy, weights)):
            da.add_unrelaxed_step(a3_mut, desc_mut)
            a3 = a3_mut
    pbs_run.relax(a3)

write('all_candidates.traj', da.get_all_relaxed_candidates())
