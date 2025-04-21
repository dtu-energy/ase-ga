from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Type, Union

import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase.neighborlist import NeighborList

__author__ = 'Stefan Bringuier <stefanbringuier@gmail.com>'
__description__ = 'LAMMPS-style native Tersoff potential for ASE'

_IMPLEMENTED_PROPERTIES = ['free_energy', 'energy', 'forces', 'stress']

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

    implemented_properties = _IMPLEMENTED_PROPERTIES

    def __init__(
        self,
        parameters: Dict[Tuple[str, str, str], TersoffParameters],
        skin: float = 0.3,
        **kwargs,
    ):
        """
        Initialize a Tersoff calculator.

        Parameters
        ----------
        parameters : dict
            Dictionary mapping element combinations to
            TersoffParameters objects.
            Format: {
                ('A', 'B', 'C'): TersoffParameters(
                    m, gamma, lambda3, c, d, h, n,
                    beta, lambda2, B, R, D, lambda1, A),
                ...
            }
            where ('A', 'B', 'C') represents the elements
            involved in the interaction.
        skin : float, default 0.3
            The skin distance for neighbor list calculations.
        **kwargs : dict
            Additional parameters to be passed to the
            ASE Calculator constructor.
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
        """
        Initialize a Tersoff calculator from a LAMMPS-style\
        Tersoff potential file.

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
        Tersoff
            Initialized Tersoff calculator with parameters from the file.
        """
        parameters = cls.read_lammps_format(potential_file)
        return cls(parameters=parameters, skin=skin, **kwargs)

    @staticmethod
    def read_lammps_format(
        potential_file: Union[str, Path],
    ) -> Dict[Tuple[str, str, str], TersoffParameters]:
        """
        Read the Tersoff potential parameters from a LAMMPS-style file.

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
        with open(potential_file, 'r') as fd:
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
    ):
        """
        Update parameters for a specific element combination.

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

    def update_nl(self, atoms) -> None:
        """
        Update the neighbor list with the parameter R+D cutoffs.

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
    ):
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
            self.update_nl(atoms)

        self.results = {}
        energy = 0.0e0
        forces = np.zeros((len(atoms), 3), dtype=np.float64)
        # Accumulated virial stress tensor
        virial = np.zeros((3, 3))

        # Duplicates atoms.get_distances() functionality, but uses
        # neighbor list's pre-computed offsets for efficiency in a
        # tight force-calculation loop rather than recompute MIC
        for i, position in enumerate(atoms.positions):
            indices, offsets = self.nl.get_neighbors(i)
            vectors = (
                atoms.positions[indices]
                + np.dot(offsets, atoms.get_cell())
                - position
            )
            distances = np.linalg.norm(vectors, axis=1)

            energy_i, force_i, stress_i = self.calc_atom_contribution(
                i, indices, distances, vectors
            )
            energy += energy_i
            forces[i] = force_i
            virial += stress_i

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces
        # Virial to stress (i.e., eV/A^3)
        stress_tensor = virial / self.atoms.get_volume()
        self.results['stress'] = stress_tensor.flat[[0, 4, 8, 5, 2, 1]]

    def calc_atom_contribution(self, i, neighbors, distances, vectors):
        """
        Calculate the energy and forces of a single atom.

        This function calculates the energy, force, and stress on atom i
        by looking at i-j pair interactions and the modification made by
        the bond order term bij with includes 3-body interaction i-j-k.

        Parameters
        ----------
        i: int
            Index of the atom
        neighbors: array_like
            Indices of the neighbor atoms
        distances: array_like
            Distances between the current atom and the neighbor atoms
        vectors: array_like
            Vectors from the current atom to the neighbor atoms

        Returns
        -------
        energy: float
            Energy contribution of the atom
        forces: array_like
            Forces on the current atom
        stress: array_like
            Stress contribution of the atom
        """
        energy = 0.0e0
        forces = np.zeros((len(neighbors), 3), dtype=np.float64)
        stress = np.zeros((3, 3), dtype=np.float64)

        type_i = self.atoms.symbols[i]

        for j, (r_ij, vec_ij) in enumerate(zip(distances, vectors)):
            type_j = self.atoms.symbols[neighbors[j]]
            key = (type_i, type_j, type_j)
            params = self.parameters[key]

            fc = self.cutoff_func(r_ij, params.R, params.D)
            if fc == 0.0:
                continue

            bij = self.calc_bond_order(j, neighbors, distances, vectors, params)

            repulsive = params.A * np.exp(-params.lambda1 * r_ij)
            attractive = -params.B * np.exp(-params.lambda2 * r_ij)

            pair_energy = fc * (repulsive + bij * attractive)
            energy += 0.5 * pair_energy

            fc_deriv = self.cutoff_func_deriv(r_ij, params.R, params.D)
            rep_deriv = -params.lambda1 * repulsive
            att_deriv = -params.lambda2 * attractive

            force_ij = -(
                (
                    fc_deriv * (repulsive + bij * attractive)
                    + fc * (rep_deriv + bij * att_deriv)
                )
                * vec_ij
                / r_ij
            )

            # Forces on neighbors j are added to i at the end
            forces[j] = force_ij

            if bij > 0.0e0:
                dbij = self.calc_bond_order_derivatives(
                    j, neighbors, distances, vectors, params
                )
                forces[j] += dbij[j] * fc * attractive

            # Virial stress
            stress += 0.5 * np.outer(vec_ij, forces[j])

        force = np.sum(forces, axis=0)
        return energy, force, stress

    def calc_bond_order(self, j, neighbors, distances, vectors, params):
        """
        Calculate bond order term considering atom i's neighbors.

        The bond order between atoms i and j is calculated as the
        sum of a few terms. The first term is the radial term, which
        is the distance between the two atoms. The second term is a
        angular term, which is the angle between the bond vector and
        the vector of atom i and its neighbors j-k. The third term is
        an exponential term, which is a function of the distance between
        the atoms and the cutoff radius.

        Parameters
        ----------
        j: int
            Index of atom
        neighbors: list of int
            List of indices of atoms that are neighbors to atom i
        distances: list of float
            List of distances between atom i and its neighbors
        vectors: list of float
            List of vectors between atom i and its neighbors
        params: dict
            Dictionary of parameters for the Tersoff potential

        Returns
        -------
        bij: float
            The bond order between atoms i and j
        """
        zeta = self._calc_zeta(j, neighbors, distances, vectors, params)
        n = params.n
        return (1.0 + params.beta**n * zeta**n) ** (-1.0 / (2.0 * n))

    def _calc_bij_d(self, zeta: float, beta: float, n: float) -> float:
        """Calculate the derivative of ``bij`` with respect to ``zeta``."""
        tmp = beta * zeta
        return (
            -0.5
            * (1.0 + tmp**n) ** (-1.0 - (1.0 / (2.0 * n)))
            * (beta * tmp ** (n - 1.0))
        )

    def _calc_zeta(self, j, neighbors, distances, vectors, params):
        """Calculate ``zeta_ij``."""
        zeta = 0.0

        for k in range(len(neighbors)):
            if k == j:
                continue

            r_ik = distances[k]
            if r_ik > params.R + params.D:
                continue

            cos_theta = np.dot(vectors[j], vectors[k]) / (
                distances[j] * distances[k]
            )
            fc_ik = self.cutoff_func(r_ik, params.R, params.D)

            g_theta = self._calc_gijk(cos_theta, params)

            # Calculate the exponential for the bond order zeta term
            # This is the term that modifies the bond order based
            # on the distance between atoms i-j and i-k. Tresholds are
            # used to prevent overflow/underflow.
            arg = (params.lambda3 * (distances[j] - r_ik)) ** params.m
            if arg > _MAX_EXP_ARG:
                ex_delr = 1.0e30
            elif arg < _MIN_EXP_ARG:
                ex_delr = 0.0
            else:
                ex_delr = np.exp(arg)

            zeta += fc_ik * g_theta * ex_delr

        return zeta

    def _calc_gijk(self, cos_theta: float, params) -> float:
        r"""Calculate the angular function ``g`` for the Tersoff potential.

        .. math::
            g(\theta) = \gamma \left( 1 + \frac{c^2}{d^2}
            - \frac{c^2}{d^2 + (h - \cos \theta)^2} \right)

        where :math:`\theta` is the angle between the bond vector
        and the vector of atom i and its neighbors j-k.
        """
        c2 = params.c * params.c
        d2 = params.d * params.d
        hcth = params.h - cos_theta
        return params.gamma * (1.0 + c2 / d2 - c2 / (d2 + hcth**2))

    def _calc_gijk_d(self, cos_theta: float, params) -> float:
        """Calculate the derivative of ``g`` with respect to ``cos_theta``."""
        c2 = params.c * params.c
        d2 = params.d * params.d
        hcth = params.h - cos_theta
        numerator = -2.0 * params.gamma * c2 * hcth
        denominator = (d2 + hcth**2) ** 2
        return numerator / denominator

    def cutoff_func(self, r: float, R: float, D: float) -> float:
        """Calculate the cutoff function."""
        if r > R + D:
            return 0.0
        if r < R - D:
            return 1.0
        return 0.5 * (1.0 - np.sin(np.pi * (r - R) / (2.0 * D)))

    def cutoff_func_deriv(self, r: float, R: float, D: float) -> float:
        """Calculate cutoff function derivative with respect to ``r``."""
        if r > R + D or r < R - D:
            return 0.0
        return -0.25 * np.pi / D * np.cos(np.pi * (r - R) / (2.0 * D))

    def calc_bond_order_derivatives(
        self, j, neighbors, distances, vectors, params
    ):
        """
        Compute the derivative of the bond order term for\
        atoms i-j.

        Parameters
        ----------
        j : int
            Index of atom j
        neighbors : array_like
            List of indices of atoms that are neighbors of atom i
        distances : array_like
            Distances between atom i and each of its neighbors
        vectors : array_like
            Vectors from atom i to each of its neighbors
        params : dict
            Tersoff parameters

        Returns
        -------
        derivatives : array_like
            Derivative of the bond order for atoms i-j

        """
        derivatives = np.zeros_like(vectors, dtype=np.float64)
        zeta = 0.0e0
        dzeta_drij = np.zeros_like(vectors[j], dtype=np.float64)

        for k in range(len(neighbors)):
            if k == j:
                continue

            r_ik = distances[k]
            vec_ik = vectors[k]
            if r_ik > params.R + params.D:
                continue

            cos_theta = np.dot(vectors[j], vec_ik) / (distances[j] * r_ik)
            fc_ik = self.cutoff_func(r_ik, params.R, params.D)
            g_theta = self._calc_gijk(cos_theta, params)

            # See comment in calc_bond_order for explanation
            arg = (params.lambda3 * (distances[j] - r_ik)) ** params.m
            if arg > _MAX_EXP_ARG:
                ex_delr = 1.0e30
            elif arg < _MIN_EXP_ARG:
                ex_delr = 0.0
            else:
                ex_delr = np.exp(arg)

            zeta += fc_ik * g_theta * ex_delr
            dzeta_drij += self._calc_zeta_d(distances[j], vec_ik, params)[1]

        bij_d = self._calc_bij_d(zeta, params.beta, params.n)

        # Derivative of the bond order (bij) w.r.t r_ij
        derivatives[j] = -1.0 * bij_d * dzeta_drij

        return derivatives

    def _calc_zeta_d(self, rij, rik, params):
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

        rij_hat = rij / np.linalg.norm(rij)
        rik_hat = rik / np.linalg.norm(rik)

        fcik = self.cutoff_func(rik, params.R, params.D)
        dfcik = self.cutoff_func_deriv(rik, params.R, params.D)

        tmp = (lam3 * (rij - rik)) ** m
        if tmp > _MAX_EXP_ARG:
            ex_delr = 1.0e30
        elif tmp < _MIN_EXP_ARG:
            ex_delr = 0.0
        else:
            ex_delr = np.exp(tmp)

        ex_delr_d = m * lam3**m * (rij - rik) ** (m - 1) * ex_delr

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
