import numpy as np

from ase.calculators.emt import EMT
from ase.ga import set_raw_score
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.population import RankFitnessPopulation
from ase.ga.slab_operators import (
    CutSpliceSlabCrossover,
    RandomCompositionMutation,
    RandomSlabPermutation,
)

# Connect to the database containing all candidates
db = DataConnection('hull.db')

# Retrieve saved parameters
pop_size = db.get_param('population_size')
refs = db.get_param('reference_energies')
metals = db.get_param('metals')
lattice_constants = db.get_param('lattice_constants')


def get_mixing_energy(atoms):
    # Set the correct cell size from the lattice constant
    new_a = get_avg_lattice_constant(atoms.get_chemical_symbols())
    # Use the orthogonal fcc cell to find the current lattice constant
    current_a = atoms.cell[0][0] / np.sqrt(2)
    atoms.set_cell(atoms.cell * new_a / current_a, scale_atoms=True)

    # Calculate the energy
    atoms.calc = EMT()
    e = atoms.get_potential_energy()

    # Subtract contributions from the pure element references
    # to get the mixing energy
    syms = atoms.get_chemical_symbols()
    for m in set(syms):
        e -= syms.count(m) * refs[m]
    return e


def get_avg_lattice_constant(syms):
    a = 0.
    for m in set(syms):
        a += syms.count(m) * lattice_constants[m]
    return a / len(syms)


def get_comp(atoms):
    return atoms.get_chemical_formula()


# Specify the number of generations this script will run
num_gens = 10

# Specify the procreation operators for the algorithm and
# how often each is picked on average
# The probability for an operator is the prepended integer divided by the sum
# of integers
oclist = [(3, CutSpliceSlabCrossover()),
          (1, RandomSlabPermutation()),
          (1, RandomCompositionMutation())
          ]
operation_selector = OperationSelector(*zip(*oclist))

# Pass parameters to the population instance
# A variable_function is required to divide candidates into groups here we use
# the chemical composition
pop = RankFitnessPopulation(data_connection=db,
                            population_size=pop_size,
                            variable_function=get_comp)

# Evaluate the starting population
# The only requirement of the evaluation is to set the raw_score
# Negative mixing energy means more stable than the pure slabs
# The optimization always progress towards larger raw score,
# so we take the negative mixing energy as the raw score
print('Evaluating initial candidates')
while db.get_number_of_unrelaxed_candidates() > 0:
    a = db.get_an_unrelaxed_candidate()
    set_raw_score(a, -get_mixing_energy(a))
    db.add_relaxed_step(a)
pop.update()

# Below is the iterative part of the algorithm
gen_num = db.get_generation_number()
for i in range(num_gens):
    print(f'Creating and evaluating generation {gen_num + i}')
    new_generation = []
    for _ in range(pop_size):
        # Select parents for a new candidate
        parents = pop.get_two_candidates()

        # Select an operator and use it
        op = operation_selector.get_operator()
        offspring, desc = op.get_new_individual(parents)
        # An operator could return None if an offspring cannot be formed
        # by the chosen parents
        if offspring is None:
            continue

        set_raw_score(offspring, -get_mixing_energy(offspring))
        new_generation.append(offspring)

    # We add a full relaxed generation at once, this is faster than adding
    # one at a time
    db.add_more_relaxed_candidates(new_generation)

    # update the population to allow new candidates to enter
    pop.update()
