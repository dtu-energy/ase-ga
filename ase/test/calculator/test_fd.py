"""Tests for `FiniteDifferenceCalculator`."""

import numpy as np
import pytest

from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.fd import FiniteDifferenceCalculator


@pytest.fixture(name='atoms')
def fixture_atoms() -> Atoms:
    """Make a fixture of atoms."""
    atoms = bulk('Cu', cubic=True)
    atoms.rattle(0.1)
    return atoms


def test_fd(atoms: Atoms) -> None:
    """Test `FiniteDifferenceCalculator`."""
    atoms.calc = EMT()
    energy_analytical = atoms.get_potential_energy()
    forces_analytical = atoms.get_forces()
    stress_analytical = atoms.get_stress()

    atoms.calc = FiniteDifferenceCalculator(EMT())
    energy_numerical = atoms.get_potential_energy()
    forces_numerical = atoms.get_forces()
    stress_numerical = atoms.get_stress()

    # check if `numerical` energy is exactly equal to `analytical`
    assert energy_numerical == energy_analytical

    # check if `numerical` forces are *not* exactly equal to `analytical`
    assert np.any(forces_numerical != forces_analytical)
    assert np.any(stress_numerical != stress_analytical)

    np.testing.assert_allclose(forces_numerical, forces_analytical)
    np.testing.assert_allclose(stress_numerical, stress_analytical)


def test_analytical_forces(atoms: Atoms) -> None:
    """Test if analytical forces are available."""
    atoms.calc = EMT()
    forces_ref = atoms.get_forces()

    atoms.calc = FiniteDifferenceCalculator(EMT(), eps_disp=None)
    forces_fdc = atoms.get_forces()

    # check if forces are *exactly* equal to `analytical`
    np.testing.assert_array_equal(forces_fdc, forces_ref)


def test_analytical_stress(atoms: Atoms) -> None:
    """Test if analytical stress is available."""
    atoms.calc = EMT()
    stress_ref = atoms.get_stress()

    atoms.calc = FiniteDifferenceCalculator(EMT(), eps_strain=None)
    stress_fdc = atoms.get_stress()

    # check if stress is *exactly* equal to `analytical`
    np.testing.assert_array_equal(stress_fdc, stress_ref)
