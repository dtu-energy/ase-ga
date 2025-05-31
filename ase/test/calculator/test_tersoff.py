"""Tests for ``Tersoff``."""

import numpy as np
import pytest

from ase import Atoms
from ase.build import bulk
from ase.calculators.calculator import PropertyNotImplementedError
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


@pytest.fixture(name='atoms_si')
def fixture_atoms_si(
    si_parameters: dict[tuple[str, str, str], TersoffParameters],
) -> Atoms:
    """Make Atoms for Si with a small displacement on the first atom."""
    atoms = bulk('Si', a=5.43, cubic=True)

    # pertubate first atom to get substantial forces
    atoms.positions[0] += [0.03, 0.02, 0.01]

    atoms.calc = Tersoff(si_parameters)

    return atoms


def test_initialize_from_params_from_dict(si_parameters):
    """Test initializing Tersoff calculator from dictionary of parameters."""
    calc = Tersoff(si_parameters)
    assert calc.parameters == si_parameters
    diamond = bulk('Si', 'diamond', a=5.43)
    diamond.calc = calc
    potential_energy = diamond.get_potential_energy()
    np.testing.assert_allclose(potential_energy, -9.260818674314585, atol=1e-8)


def test_set_parameters(
    si_parameters: dict[tuple[str, str, str], TersoffParameters],
) -> None:
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


def test_isolated_atom(si_parameters: dict) -> None:
    """Test if an isolated atom can be computed correctly."""
    atoms = Atoms('Si')
    atoms.calc = Tersoff(si_parameters)
    energy = atoms.get_potential_energy()
    energies = atoms.get_potential_energies()
    forces = atoms.get_forces()
    np.testing.assert_almost_equal(energy, 0.0)
    np.testing.assert_allclose(energies, [0.0], rtol=1e-5)
    np.testing.assert_allclose(forces, [[0.0] * 3], rtol=1e-5)
    with pytest.raises(PropertyNotImplementedError):
        atoms.get_stress()


def test_unary(atoms_si: Atoms) -> None:
    """Test if energy, forces, and stress of a unary system agree with LAMMPS.

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
    energy_ref = -37.03237572778589
    energies_ref = [
        -4.62508202,
        -4.62242901,
        -4.63032346,
        -4.63028909,
        -4.63037555,
        -4.63147495,
        -4.63040683,
        -4.63199482,
    ]
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

    energy = atoms_si.get_potential_energy()
    energies = atoms_si.get_potential_energies()
    forces = atoms_si.get_forces()
    stress = atoms_si.get_stress()
    np.testing.assert_almost_equal(energy, energy_ref)
    np.testing.assert_allclose(energies, energies_ref, rtol=1e-5)
    np.testing.assert_allclose(forces, forces_ref, rtol=1e-5)
    np.testing.assert_allclose(stress, stress_ref, rtol=1e-5)


def test_binary(datadir) -> None:
    """Test if energy, forces, and stress of a binary system agree with LAMMPS.

    The reference values are obtained in the following way.

    >>> from ase.calculators.lammpslib import LAMMPSlib
    >>>
    >>> atoms = bulk('Si', a=5.43, cubic=True)
    >>> atoms.symbols[1] = 'C'
    >>> atoms.symbols[2] = 'C'
    >>> atoms.positions[0] += [0.03, 0.02, 0.01]
    >>> lmpcmds = ['pair_style tersoff', 'pair_coeff * * SiC.tersoff Si C']
    >>> atoms.calc = LAMMPSlib(lmpcmds=lmpcmds)
    >>> energy = atoms.get_potential_energy()
    >>> energies = atoms.get_potential_energies()
    >>> forces = atoms.get_forces()
    >>> stress = atoms.get_stress()

    """
    atoms = bulk('Si', a=5.43, cubic=True)
    atoms.symbols[1] = 'C'
    atoms.symbols[2] = 'C'

    # pertubate first atom to get substantial forces
    atoms.positions[0] += [0.03, 0.02, 0.01]

    potential_file = datadir / 'tersoff' / 'SiC.tersoff'
    atoms.calc = Tersoff.from_lammps(potential_file)

    energy_ref = -28.780184609451915
    energies_ref = [
        -4.33637575,
        -2.02218449,
        -1.80044260,
        -4.12192108,
        -4.12650203,
        -4.12473794,
        -4.12677193,
        -4.12124880,
    ]
    forces_ref = [
        [+6.40479511, +6.64830387, +6.83733140],
        [+6.93259841, -7.29932178, -7.32986722],
        [-7.09646214, +7.17384006, +7.15478087],
        [-6.63906798, -6.62943844, -6.64034900],
        [-6.48323569, +6.55237593, -6.52463929],
        [+6.64668748, +6.56474714, -6.55390547],
        [-6.47026652, -6.50615747, +6.53977908],
        [+6.70495132, -6.50434930, +6.51686963],
    ]
    stress_ref = [
        +0.35188635,
        +0.35366730,
        +0.35444532,
        -0.11629806,
        +0.11665436,
        +0.11668285,
    ]

    energy = atoms.get_potential_energy()
    energies = atoms.get_potential_energies()
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    np.testing.assert_almost_equal(energy, energy_ref)
    np.testing.assert_allclose(energies, energies_ref, rtol=1e-5)
    np.testing.assert_allclose(forces, forces_ref, rtol=1e-5)
    np.testing.assert_allclose(stress, stress_ref, rtol=1e-5)


def test_forces_and_stress(atoms_si: Atoms) -> None:
    """Test if analytical forces and stress agree with numerical ones."""
    forces = atoms_si.get_forces()
    numerical_forces = calculate_numerical_forces(atoms_si, eps=1e-5)
    np.testing.assert_allclose(forces, numerical_forces, atol=1e-5)

    stress = atoms_si.get_stress()
    numerical_stress = calculate_numerical_stress(atoms_si, eps=1e-5)
    np.testing.assert_allclose(stress, numerical_stress, atol=1e-5)
