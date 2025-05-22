from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Type, Union

import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase.neighborlist import NeighborList
from ase.stress import full_3x3_to_voigt_6_stress

__author__ = 'Stefan Bringuier <stefanbringuier@gmail.com>'
__description__ = 'LAMMPS-style native Tersoff potential for ASE'

# Maximum/minimum exponents for numerical stability
# in bond order calculation
_MAX_EXP_ARG = 69.0776e0
_MIN_EXP_ARG = -69.0776e0


@dataclass
class TersoffParameters:
    """Parameters for 3 element Tersoff potential interaction.

    Can be instantiated with either positional or keyword arguments:
        TersoffParameters(1.0, 2.0, ...) or
        TersoffParameters(m=1.0, gamma=2.0, ...)
    """

    m: float
    gamma: float
    lambda3: float
    c: float
    d: float
    h: float
    n: float
    beta: float
    lambda2: float
    B: float
    R: float
    D: float
    lambda1: float
    A: float

    @classmethod
    def from_list(cls, params: List[float]) -> 'TersoffParameters':
        """Create TersoffParameters from a list of 14 parameter values."""
        if len(params) != 14:
            raise ValueError(f'Expected 14 parameters, got {len(params)}')
        return cls(*map(float, params))


class Tersoff(Calculator):
    """ASE Calculator for Tersoff interatomic potential.

    .. versionadded:: 3.25.0
    """

    implemented_properties = [
        'free_energy',
        'energy',
        'energies',
        'forces',
        'stress',
    ]

    def __init__(
        self,
        parameters: Dict[Tuple[str, str, str], TersoffParameters],
        skin: float = 0.3,
        **kwargs,
    ) -> None:
        """
        Parameters
        ----------
        parameters : dict
            Mapping element combinations to TersoffParameters objects::

                {
                    ('A', 'B', 'C'): TersoffParameters(
                        m, gamma, lambda3, c, d, h, n,
                        beta, lambda2, B, R, D, lambda1, A),
                    ...
                }

            where ('A', 'B', 'C') are the elements involved in the interaction.
        skin : float, default 0.3
            The skin distance for neighbor list calculations.
        **kwargs : dict
            Additional parameters to be passed to
            :class:`~ase.calculators.Calculator`.

        """
        Calculator.__init__(self, **kwargs)
        self.cutoff_skin = skin
        self.parameters = parameters

    @classmethod
    def from_lammps(
        cls: Type['Tersoff'],
        potential_file: Union[str, Path],
        skin: float = 0.3,
        **kwargs,
    ) -> 'Tersoff':
        """Make :class:`Tersoff` from a LAMMPS-style Tersoff potential file.

        Parameters
        ----------
        potential_file : str or Path
            The path to a LAMMPS-style Tersoff potential file.
        skin : float, default 0.3
            The skin distance for neighbor list calculations.
        **kwargs : dict
            Additional parameters to be passed to the
            ASE Calculator constructor.

        Returns
        -------
        :class:`Tersoff`
            Initialized Tersoff calculator with parameters from the file.

        """
        parameters = cls.read_lammps_format(potential_file)
        return cls(parameters=parameters, skin=skin, **kwargs)

    @staticmethod
    def read_lammps_format(
        potential_file: Union[str, Path],
    ) -> Dict[Tuple[str, str, str], TersoffParameters]:
        """Read the Tersoff potential parameters from a LAMMPS-style file.

        Parameters
        ----------
        potential_file : str or Path
            Path to the LAMMPS-style Tersoff potential file

        Returns
        -------
        dict
            Dictionary mapping element combinations to TersoffParameters objects

        """
        block_size = 17
        with Path(potential_file).open('r', encoding='utf-8') as fd:
            content = (
                ''.join(
                    [line for line in fd if not line.strip().startswith('#')]
                )
                .replace('\n', ' ')
                .split()
            )

        if len(content) % block_size != 0:
            raise ValueError(
                'The potential file does not have the correct LAMMPS format.'
            )

        parameters: Dict[Tuple[str, str, str], TersoffParameters] = {}
        for i in range(0, len(content), block_size):
            block = content[i : i + block_size]
            e1, e2, e3 = block[0], block[1], block[2]
            current_elements = (e1, e2, e3)
            params = map(float, block[3:])
            parameters[current_elements] = TersoffParameters(*params)

        return parameters

    def set_parameters(
        self,
        key: Tuple[str, str, str],
        params: TersoffParameters = None,
        **kwargs,
    ) -> None:
        """Update parameters for a specific element combination.

        Parameters
        ----------
        key: Tuple[str, str, str]
            The element combination key of the parameters to be updated
        params: TersoffParameters, optional
            A TersoffParameters instance to completely replace the parameters
        **kwargs:
            Individual parameter values to update, e.g. R=2.9

        """
        if key not in self.parameters:
            raise KeyError(f"Key '{key}' not found in parameters.")

        if params is not None:
            if kwargs:
                raise ValueError('Cannot provide both params and kwargs.')
            self.parameters[key] = params
        else:
            for name, value in kwargs.items():
                if not hasattr(self.parameters[key], name):
                    raise ValueError(f'Invalid parameter name: {name}')
                setattr(self.parameters[key], name, value)

    def _update_nl(self, atoms) -> None:
        """Update the neighbor list with the parameter R+D cutoffs.

        Parameters
        ----------
        atoms: ase.Atoms
            The atoms to calculate the neighbor list for.

        Notes
        -----
        The cutoffs are determined by the parameters of the Tersoff potential.
        Each atom's cutoff is based on the R+D values from the parameter set
        where that atom's element appears first in the key tuple.

        """
        # Get cutoff for each atom based on its element type
        cutoffs = []

        for symbol in atoms.symbols:
            # Find first parameter set, element is the first slot
            param_key = next(
                key for key in self.parameters.keys() if key[0] == symbol
            )
            params = self.parameters[param_key]
            cutoffs.append(params.R + params.D)

        self.nl = NeighborList(
            cutoffs,
            skin=self.cutoff_skin,
            self_interaction=False,
            bothways=True,
        )

        self.nl.update(atoms)

    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ) -> None:
        """Calculate energy, forces, and stress.

        Notes
        -----
        The force and stress are calculated regardless if they are
        requested, despite some additional overhead cost,
        therefore they are always stored in the results dict.

        """
        Calculator.calculate(self, atoms, properties, system_changes)

        # Rebuild neighbor list when any relevant system changes occur
        checks = {'positions', 'numbers', 'cell', 'pbc'}
        if any(change in checks for change in system_changes) or not hasattr(
            self, 'nl'
        ):
            self._update_nl(atoms)

        self.results = {}
        energies = np.zeros(len(atoms))
        forces = np.zeros((len(atoms), 3))
        virial = np.zeros((3, 3))

        # Duplicates atoms.get_distances() functionality, but uses
        # neighbor list's pre-computed offsets for efficiency in a
        # tight force-calculation loop rather than recompute MIC
        for i in range(len(atoms)):
            self._calc_atom_contribution(i, energies, forces, virial)

        self.results['energies'] = energies
        self.results['energy'] = self.results['free_energy'] = energies.sum()
        self.results['forces'] = forces
        # Virial to stress (i.e., eV/A^3)
        if self.atoms.cell.rank == 3:
            stress = virial / self.atoms.get_volume()
            self.results['stress'] = full_3x3_to_voigt_6_stress(stress)

    def _calc_atom_contribution(
        self,
        idx_i: int,
        energies: np.ndarray,
        forces: np.ndarray,
        virial: np.ndarray,
    ) -> None:
        """Calculate the contributions of a single atom to the properties.

        This function calculates the energy, force, and stress on atom i
        by looking at i-j pair interactions and the modification made by
        the bond order term bij with includes 3-body interaction i-j-k.

        Parameters
        ----------
        idx_i: int
            Index of atom i
        energies: array_like
            Site energies to be updated.
        forces: array_like
            Forces to be updated.
        virial: array_like
            Virial tensor to be updated.

        """
        indices, offsets = self.nl.get_neighbors(idx_i)
        vectors = self.atoms.positions[indices]
        vectors += offsets @ self.atoms.cell
        vectors -= self.atoms.positions[idx_i]
        distances = np.sqrt(np.add.reduce(vectors**2, axis=1))

        type_i = self.atoms.symbols[idx_i]
        for j, (idx_j, abs_rij, rij) in enumerate(
            zip(indices, distances, vectors)
        ):
            type_j = self.atoms.symbols[idx_j]
            key = (type_i, type_j, type_j)
            params = self.parameters[key]

            rij_hat = rij / abs_rij

            fc = self._calc_fc(abs_rij, params.R, params.D)
            if fc == 0.0:
                continue

            zeta = self._calc_zeta(type_i, j, indices, distances, vectors)
            bij = self._calc_bij(zeta, params.beta, params.n)
            bij_d = self._calc_bij_d(zeta, params.beta, params.n)

            repulsive = params.A * np.exp(-params.lambda1 * abs_rij)
            attractive = -params.B * np.exp(-params.lambda2 * abs_rij)

            # distribute the pair energy evenly to be consistent with LAMMPS
            energies[idx_i] += 0.25 * fc * (repulsive + bij * attractive)
            energies[idx_j] += 0.25 * fc * (repulsive + bij * attractive)

            dfc = self._calc_fc_d(abs_rij, params.R, params.D)
            rep_deriv = -params.lambda1 * repulsive
            att_deriv = -params.lambda2 * attractive

            tmp = dfc * (repulsive + bij * attractive)
            tmp += fc * (rep_deriv + bij * att_deriv)

            # derivative with respect to the position of atom j
            grad = 0.5 * tmp * rij_hat

            forces[idx_i] += grad
            forces[idx_j] -= grad

            virial += np.outer(grad, rij)

            for k, idx_k in enumerate(indices):
                if k == j:
                    continue

                type_k = self.atoms.symbols[idx_k]
                key = (type_i, type_j, type_k)
                params = self.parameters[key]

                if distances[k] > params.R + params.D:
                    continue

                rik = vectors[k]

                dztdri, dztdrj, dztdrk = self._calc_zeta_d(rij, rik, params)

                gradi = 0.5 * fc * bij_d * dztdri * attractive
                gradj = 0.5 * fc * bij_d * dztdrj * attractive
                gradk = 0.5 * fc * bij_d * dztdrk * attractive

                forces[idx_i] -= gradi
                forces[idx_j] -= gradj
                forces[idx_k] -= gradk

                virial += np.outer(gradj, rij)
                virial += np.outer(gradk, rik)

    def _calc_bij(self, zeta: float, beta: float, n: float) -> float:
        """Calculate the bond order ``bij`` between atoms ``i`` and ``j``."""
        tmp = beta * zeta
        return (1.0 + tmp**n) ** (-1.0 / (2.0 * n))

    def _calc_bij_d(self, zeta: float, beta: float, n: float) -> float:
        """Calculate the derivative of ``bij`` with respect to ``zeta``."""
        tmp = beta * zeta
        return (
            -0.5
            * (1.0 + tmp**n) ** (-1.0 - (1.0 / (2.0 * n)))
            * (beta * tmp ** (n - 1.0))
        )

    def _calc_zeta(
        self,
        type_i: str,
        j: int,
        neighbors: np.ndarray,
        distances: np.ndarray,
        vectors: np.ndarray,
    ) -> float:
        """Calculate ``zeta_ij``."""
        idx_j = neighbors[j]
        type_j = self.atoms.symbols[idx_j]
        abs_rij = distances[j]

        zeta = 0.0

        for k, idx_k in enumerate(neighbors):
            if k == j:
                continue

            type_k = self.atoms.symbols[idx_k]
            key = (type_i, type_j, type_k)
            params = self.parameters[key]

            abs_rik = distances[k]
            if abs_rik > params.R + params.D:
                continue

            costheta = np.dot(vectors[j], vectors[k]) / (abs_rij * abs_rik)
            fc_ik = self._calc_fc(abs_rik, params.R, params.D)

            g_theta = self._calc_gijk(costheta, params)

            # Calculate the exponential for the bond order zeta term
            # This is the term that modifies the bond order based
            # on the distance between atoms i-j and i-k. Tresholds are
            # used to prevent overflow/underflow.
            arg = (params.lambda3 * (abs_rij - abs_rik)) ** params.m
            if arg > _MAX_EXP_ARG:
                ex_delr = 1.0e30
            elif arg < _MIN_EXP_ARG:
                ex_delr = 0.0
            else:
                ex_delr = np.exp(arg)

            zeta += fc_ik * g_theta * ex_delr

        return zeta

    def _calc_gijk(self, costheta: float, params: TersoffParameters) -> float:
        r"""Calculate the angular function ``g`` for the Tersoff potential.

        .. math::
            g(\theta) = \gamma \left( 1 + \frac{c^2}{d^2}
            - \frac{c^2}{d^2 + (h - \cos \theta)^2} \right)

        where :math:`\theta` is the angle between the bond vector
        and the vector of atom i and its neighbors j-k.
        """
        c2 = params.c * params.c
        d2 = params.d * params.d
        hcth = params.h - costheta
        return params.gamma * (1.0 + c2 / d2 - c2 / (d2 + hcth**2))

    def _calc_gijk_d(self, costheta: float, params: TersoffParameters) -> float:
        """Calculate the derivative of ``g`` with respect to ``costheta``."""
        c2 = params.c * params.c
        d2 = params.d * params.d
        hcth = params.h - costheta
        numerator = -2.0 * params.gamma * c2 * hcth
        denominator = (d2 + hcth**2) ** 2
        return numerator / denominator

    def _calc_fc(self, r: np.floating, R: float, D: float) -> float:
        """Calculate the cutoff function."""
        if r > R + D:
            return 0.0
        if r < R - D:
            return 1.0
        return 0.5 * (1.0 - np.sin(np.pi * (r - R) / (2.0 * D)))

    def _calc_fc_d(self, r: np.floating, R: float, D: float) -> float:
        """Calculate cutoff function derivative with respect to ``r``."""
        if r > R + D or r < R - D:
            return 0.0
        return -0.25 * np.pi / D * np.cos(np.pi * (r - R) / (2.0 * D))

    def _calc_zeta_d(
        self,
        rij: np.ndarray,
        rik: np.ndarray,
        params: TersoffParameters,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate the derivatives of ``zeta``.

        Returns
        -------
        dri : ndarray of shape (3,), dtype float
            Derivative with respect to the position of atom ``i``.
        drj : ndarray of shape (3,), dtype float
            Derivative with respect to the position of atom ``j``.
        drk : ndarray of shape (3,), dtype float
            Derivative with respect to the position of atom ``k``.

        """
        lam3 = params.lambda3
        m = params.m

        abs_rij = np.linalg.norm(rij)
        abs_rik = np.linalg.norm(rik)

        rij_hat = rij / abs_rij
        rik_hat = rik / abs_rik

        fcik = self._calc_fc(abs_rik, params.R, params.D)
        dfcik = self._calc_fc_d(abs_rik, params.R, params.D)

        tmp = (lam3 * (abs_rij - abs_rik)) ** m
        if tmp > _MAX_EXP_ARG:
            ex_delr = 1.0e30
        elif tmp < _MIN_EXP_ARG:
            ex_delr = 0.0
        else:
            ex_delr = np.exp(tmp)

        ex_delr_d = m * lam3**m * (abs_rij - abs_rik) ** (m - 1) * ex_delr

        costheta = rij_hat @ rik_hat
        gijk = self._calc_gijk(costheta, params)
        gijk_d = self._calc_gijk_d(costheta, params)

        dcosdri, dcosdrj, dcosdrk = self._calc_costheta_d(rij, rik)

        dri = -dfcik * gijk * ex_delr * rik_hat
        dri += fcik * gijk_d * ex_delr * dcosdri
        dri += fcik * gijk * ex_delr_d * rik_hat
        dri -= fcik * gijk * ex_delr_d * rij_hat

        drj = fcik * gijk_d * ex_delr * dcosdrj
        drj += fcik * gijk * ex_delr_d * rij_hat

        drk = dfcik * gijk * ex_delr * rik_hat
        drk += fcik * gijk_d * ex_delr * dcosdrk
        drk -= fcik * gijk * ex_delr_d * rik_hat

        return dri, drj, drk

    def _calc_costheta_d(
        self,
        rij: np.ndarray,
        rik: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""Calculate the derivatives of ``costheta``.

        If

        .. math::
            \cos \theta = \frac{\mathbf{u} \cdot \mathbf{v}}{u v}

        Then

        .. math::
            \frac{\partial \cos \theta}{\partial \mathbf{u}}
            = \frac{\mathbf{v}}{u v}
            - \frac{\mathbf{u} \cdot \mathbf{v}}{v} \cdot \frac{\mathbf{u}}{u^3}
            = \frac{\mathbf{v}}{u v} - \frac{\cos \theta}{u^2} \mathbf{u}

        Parameters
        ----------
        rij : ndarray of shape (3,), dtype float
            Vector from atoms ``i`` to ``j``.
        rik : ndarray of shape (3,), dtype float
            Vector from atoms ``i`` to ``k``.

        Returns
        -------
        dri : ndarray of shape (3,), dtype float
            Derivative with respect to the position of atom ``i``.
        drj : ndarray of shape (3,), dtype float
            Derivative with respect to the position of atom ``j``.
        drk : ndarray of shape (3,), dtype float
            Derivative with respect to the position of atom ``k``.

        """
        abs_rij = np.linalg.norm(rij)
        abs_rik = np.linalg.norm(rik)
        costheta = (rij @ rik) / (abs_rij * abs_rik)
        drj = (rik / abs_rik - costheta * rij / abs_rij) / abs_rij
        drk = (rij / abs_rij - costheta * rik / abs_rik) / abs_rik
        dri = -(drj + drk)
        return dri, drj, drk
