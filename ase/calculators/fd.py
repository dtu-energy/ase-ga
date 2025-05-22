"""Module for `FiniteDifferenceCalculator`."""

from collections.abc import Iterable
from functools import partial
from typing import Optional

import numpy as np

from ase import Atoms
from ase.calculators.calculator import BaseCalculator, all_properties


class FiniteDifferenceCalculator(BaseCalculator):
    """Wrapper calculator using the finite-difference method.

    The forces and the stress are computed using the finite-difference method.

    .. versionadded:: 3.24.0
    """

    implemented_properties = all_properties

    def __init__(
        self,
        calc: BaseCalculator,
        eps_disp: Optional[float] = 1e-6,
        eps_strain: Optional[float] = 1e-6,
        *,
        force_consistent: bool = True,
    ) -> None:
        """

        Parameters
        ----------
        calc : :class:`~ase.calculators.calculator.BaseCalculator`
            ASE Calculator object to be wrapped.
        eps_disp : Optional[float], default 1e-6
            Displacement used for computing forces.
            If :py:obj:`None`, analytical forces are computed.
        eps_strain : Optional[float], default 1e-6
            Strain used for computing stress.
            If :py:obj:`None`, analytical stress is computed.
        force_consistent : bool, default :py:obj:`True`
            If :py:obj:`True`, the energies consistent with the forces are used
            for finite-difference calculations.

        """
        super().__init__()
        self.calc = calc
        self.eps_disp = eps_disp
        self.eps_strain = eps_strain
        self.force_consistent = force_consistent

    def calculate(self, atoms: Atoms, properties, system_changes) -> None:
        atoms = atoms.copy()  # copy to not mess up original `atoms`
        atoms.calc = self.calc
        self.results = {}
        self.results['energy'] = self.calc.get_potential_energy(atoms)
        for key in ['free_energy']:
            if key in self.calc.results:
                self.results[key] = self.calc.results[key]
        if self.eps_disp is None:
            self.results['forces'] = self.calc.get_forces(atoms)
        else:
            self.results['forces'] = calculate_numerical_forces(
                atoms,
                eps=self.eps_disp,
                force_consistent=self.force_consistent,
            )
        if self.eps_strain is None:
            self.results['stress'] = self.calc.get_stress(atoms)
        else:
            self.results['stress'] = calculate_numerical_stress(
                atoms,
                eps=self.eps_strain,
                force_consistent=self.force_consistent,
            )


def _numeric_force(
    atoms: Atoms,
    iatom: int,
    icart: int,
    eps: float = 1e-6,
    *,
    force_consistent: bool = False,
) -> float:
    """Calculate numerical force on a specific atom along a specific direction.

    Parameters
    ----------
    atoms : :class:`~ase.Atoms`
        ASE :class:`~ase.Atoms` object.
    iatom : int
        Index of atoms.
    icart : {0, 1, 2}
        Index of Cartesian component.
    eps : float, default 1e-6
        Displacement.
    force_consistent : bool, default :py:obj:`False`
        If :py:obj:`True`, the energies consistent with the forces are used for
        finite-difference calculations.

    """
    p0 = atoms.get_positions()
    p = p0.copy()
    p[iatom, icart] = p0[iatom, icart] + eps
    atoms.set_positions(p, apply_constraint=False)
    eplus = atoms.get_potential_energy(force_consistent=force_consistent)
    p[iatom, icart] = p0[iatom, icart] - eps
    atoms.set_positions(p, apply_constraint=False)
    eminus = atoms.get_potential_energy(force_consistent=force_consistent)
    atoms.set_positions(p0, apply_constraint=False)
    return (eminus - eplus) / (2 * eps)


def calculate_numerical_forces(
    atoms: Atoms,
    eps: float = 1e-6,
    iatoms: Optional[Iterable[int]] = None,
    icarts: Optional[Iterable[int]] = None,
    *,
    force_consistent: bool = False,
) -> np.ndarray:
    """Calculate forces numerically based on the finite-difference method.

    Parameters
    ----------
    atoms : :class:`~ase.Atoms`
        ASE :class:`~ase.Atoms` object.
    eps : float, default 1e-6
        Displacement.
    iatoms : Optional[Iterable[int]]
        Indices of atoms for which forces are computed.
        By default, all atoms are considered.
    icarts : Optional[Iterable[int]]
        Indices of Cartesian coordinates for which forces are computed.
        By default, all three coordinates are considered.
    force_consistent : bool, default :py:obj:`False`
        If :py:obj:`True`, the energies consistent with the forces are used for
        finite-difference calculations.

    Returns
    -------
    forces : np.ndarray
        Forces computed numerically based on the finite-difference method.

    """
    if iatoms is None:
        iatoms = range(len(atoms))
    if icarts is None:
        icarts = [0, 1, 2]
    f = partial(_numeric_force, eps=eps, force_consistent=force_consistent)
    forces = [[f(atoms, iatom, icart) for icart in icarts] for iatom in iatoms]
    return np.array(forces)


def calculate_numerical_stress(
    atoms: Atoms,
    eps: float = 1e-6,
    voigt: bool = True,
    *,
    force_consistent: bool = True,
) -> np.ndarray:
    """Calculate stress numerically based on the finite-difference method.

    Parameters
    ----------
    atoms : :class:`~ase.Atoms`
        ASE :class:`~ase.Atoms` object.
    eps : float, default 1e-6
        Strain in the Voigt notation.
    voigt : bool, default :py:obj:`True`
        If :py:obj:`True`, the stress is returned in the Voigt notation.
    force_consistent : bool, default :py:obj:`True`
        If :py:obj:`True`, the energies consistent with the forces are used for
        finite-difference calculations.

    Returns
    -------
    stress : np.ndarray
        Stress computed numerically based on the finite-difference method.

    """
    stress = np.zeros((3, 3), dtype=float)

    cell = atoms.cell.copy()
    volume = atoms.get_volume()
    for i in range(3):
        x = np.eye(3)
        x[i, i] = 1.0 + eps
        atoms.set_cell(cell @ x, scale_atoms=True)
        eplus = atoms.get_potential_energy(force_consistent=force_consistent)

        x[i, i] = 1.0 - eps
        atoms.set_cell(cell @ x, scale_atoms=True)
        eminus = atoms.get_potential_energy(force_consistent=force_consistent)

        stress[i, i] = (eplus - eminus) / (2 * eps * volume)
        x[i, i] = 1.0

        j = i - 2
        x[i, j] = x[j, i] = +0.5 * eps
        atoms.set_cell(cell @ x, scale_atoms=True)
        eplus = atoms.get_potential_energy(force_consistent=force_consistent)

        x[i, j] = x[j, i] = -0.5 * eps
        atoms.set_cell(cell @ x, scale_atoms=True)
        eminus = atoms.get_potential_energy(force_consistent=force_consistent)

        stress[i, j] = stress[j, i] = (eplus - eminus) / (2 * eps * volume)

    atoms.set_cell(cell, scale_atoms=True)

    return stress.flat[[0, 4, 8, 5, 2, 1]] if voigt else stress
