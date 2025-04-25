"""Tests for ``Tersoff``."""

import numpy as np
import pytest

from ase.build import bulk
from ase.calculators.fd import (
    calculate_numerical_forces,
    calculate_numerical_stress,
)
from ase.calculators.tersoff import Tersoff, TersoffParameters


@pytest.fixture
def si_parameters():
    """Fixture providing the Silicon parameters.

    Parameters taken from: Tersoff, Phys Rev B, 37, 6991 (1988)
    """
    return {
        ('Si', 'Si', 'Si'): TersoffParameters(
            A=3264.7,
            B=95.373,
            lambda1=3.2394,
            lambda2=1.3258,
            lambda3=1.3258,
            beta=0.33675,
            gamma=1.00,
            m=3.00,
            n=22.956,
            c=4.8381,
            d=2.0417,
            h=0.0000,
            R=3.00,
            D=0.20,
        )
    }


def test_initialize_from_params_from_dict(si_parameters):
    """Test initializing Tersoff calculator from dictionary of parameters."""
    calc = Tersoff(si_parameters)
    assert calc.parameters == si_parameters
    diamond = bulk('Si', 'diamond', a=5.43)
    diamond.calc = calc
    potential_energy = diamond.get_potential_energy()
    np.testing.assert_allclose(potential_energy, -9.260818674314585, atol=1e-8)


def test_set_parameters(si_parameters: dict[tuple, TersoffParameters]) -> None:
    """Test updating parameters of the Tersoff calculator."""
    calc = Tersoff(si_parameters)
    key = ('Si', 'Si', 'Si')

    calc.set_parameters(key, m=2.0)
    assert calc.parameters[key].m == 2.0

    calc.set_parameters(key, R=2.90, D=0.25)
    assert calc.parameters[key].R == 2.90
    assert calc.parameters[key].D == 0.25

    new_params = TersoffParameters(
        m=si_parameters[key].m,
        gamma=si_parameters[key].gamma,
        lambda3=si_parameters[key].lambda3,
        c=si_parameters[key].c,
        d=si_parameters[key].d,
        h=si_parameters[key].h,
        n=si_parameters[key].n,
        beta=si_parameters[key].beta,
        lambda2=si_parameters[key].lambda2,
        B=si_parameters[key].B,
        R=3.00,  # Reset cutoff radius
        D=si_parameters[key].D,
        lambda1=si_parameters[key].lambda1,
        A=si_parameters[key].A,
    )
    calc.set_parameters(key, params=new_params)
    assert calc.parameters[key] == new_params


def test_properties(si_parameters: dict) -> None:
    """Test if energy, forces, and stress agree with LAMMPS.

    The reference values are obtained in the following way.

    >>> from ase.calculators.lammpslib import LAMMPSlib
    >>>
    >>> atoms = bulk('Si', a=5.43, cubic=True)
    >>> atoms.positions[0] += [0.03, 0.02, 0.01]
    >>> lmpcmds = ['pair_style tersoff', 'pair_coeff * * Si.tersoff Si']
    >>> atoms.calc = LAMMPSlib(lmpcmds=lmpcmds)
    >>> energy = atoms.get_potential_energy()
    >>> energies = atoms.get_potential_energies()
    >>> forces = atoms.get_forces()
    >>> stress = atoms.get_stress()

    """
    atoms = bulk('Si', a=5.43, cubic=True)

    # pertubate first atom to get substantial forces
    atoms.positions[0] += [0.03, 0.02, 0.01]

    atoms.calc = Tersoff(si_parameters)

    energy_ref = -37.03237572778589
    forces_ref = [
        [-4.63805736e-01, -3.17112011e-01, -1.79345801e-01],
        [+2.34142607e-01, +2.29060580e-01, +2.24142706e-01],
        [-2.79544489e-02, +1.31289732e-03, +3.99485914e-04],
        [+1.85144670e-02, +1.48017753e-02, +8.47421196e-03],
        [+2.06558877e-03, -1.86613107e-02, +3.98039278e-04],
        [+8.68756690e-02, -5.15405628e-02, +7.32472691e-02],
        [+2.06388309e-03, +1.30960793e-03, -9.30764103e-03],
        [+1.48097970e-01, +1.40829025e-01, -1.18008270e-01],
    ]
    stress_ref = [
        -0.00048610,
        -0.00056779,
        -0.00061684,
        -0.00342602,
        -0.00231541,
        -0.00124569,
    ]

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    np.testing.assert_almost_equal(energy, energy_ref)
    np.testing.assert_allclose(forces, forces_ref, rtol=1e-5)
    np.testing.assert_allclose(stress, stress_ref, rtol=1e-5)


def test_forces_and_stress(si_parameters: dict) -> None:
    """Test if analytical forces and stress agree with numerical ones."""
    atoms = bulk('Si', a=5.43, cubic=True)

    # pertubate first atom to get substantial forces
    atoms.positions[0] += [0.03, 0.02, 0.01]

    atoms.calc = Tersoff(si_parameters)

    forces = atoms.get_forces()
    numerical_forces = calculate_numerical_forces(atoms, eps=1e-5)
    np.testing.assert_allclose(forces, numerical_forces, atol=1e-5)

    stress = atoms.get_stress()
    numerical_stress = calculate_numerical_stress(atoms, eps=1e-5)
    np.testing.assert_allclose(stress, numerical_stress, atol=1e-5)
