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
