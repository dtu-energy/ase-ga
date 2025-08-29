import numpy as np

from ase.calculators.emt import EMT
from ase.db import connect
from ase.eos import EquationOfState
from ase.lattice.cubic import FaceCenteredCubic


def relax(input_atoms, ref_db):
    atoms_string = input_atoms.get_chemical_symbols()

    # Open connection to the database with reference data
    db = connect(ref_db)

    # Load our model structure which is just FCC
    atoms = FaceCenteredCubic('X', latticeconstant=1.0)
    atoms.set_chemical_symbols(atoms_string)

    # Compute the average lattice constant of the metals in this individual
    # and the sum of energies of the constituent metals in the fcc lattice
    # we will need this for calculating the heat of formation
    a = 0
    ei = 0
    for m in set(atoms_string):
        dct = db.get(metal=m)
        count = atoms_string.count(m)
        a += count * dct.latticeconstant
        ei += count * dct.energy_per_atom
    a /= len(atoms_string)
    atoms.set_cell([a, a, a], scale_atoms=True)

    # Since calculations are extremely fast with EMT we can also do a volume
    # relaxation
    atoms.calc = EMT()
    eps = 0.05
    volumes = (a * np.linspace(1 - eps, 1 + eps, 9)) ** 3
    energies = []
    for v in volumes:
        atoms.set_cell([v ** (1.0 / 3)] * 3, scale_atoms=True)
        energies.append(atoms.get_potential_energy())

    eos = EquationOfState(volumes, energies)
    v1, ef, _B = eos.fit()
    latticeconstant = v1 ** (1.0 / 3)

    # Calculate the heat of formation by subtracting ef with ei
    hof = (ef - ei) / len(atoms)

    # Place the calculated parameters in the info dictionary of the
    # input_atoms object
    input_atoms.info['key_value_pairs']['hof'] = hof

    # Raw score must always be set
    # Use one of the following two; they are equivalent
    input_atoms.info['key_value_pairs']['raw_score'] = -hof
    # set_raw_score(input_atoms, -hof)

    input_atoms.info['key_value_pairs']['latticeconstant'] = latticeconstant

    # Setting the atoms_string directly for easier analysis
    atoms_string = ''.join(input_atoms.get_chemical_symbols())
    input_atoms.info['key_value_pairs']['atoms_string'] = atoms_string
