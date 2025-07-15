from abc import ABC, abstractmethod

import numpy as np

# Due to the high prevalence of cyclic imports surrounding ase.optimize,
# we define the Optimizable ABC here in utils.
# Can we find a better way?


class Optimizable(ABC):
    @abstractmethod
    def ndofs(self) -> int:
        """Return number of degrees of freedom."""

    @abstractmethod
    def get_x(self) -> np.ndarray:
        """Return current coordinates as a flat ndarray."""

    @abstractmethod
    def set_x(self, x: np.ndarray) -> None:
        """Set flat ndarray as current coordinates."""

    @abstractmethod
    def get_gradient(self) -> np.ndarray:
        """Return gradient at current coordinates as flat ndarray.

        NOTE: Currently this is still the (flat) "forces" i.e.
        the negative gradient.  This must be fixed before the optimizable
        API is done."""
        # Callers who want Nx3 will do ".get_gradient().reshape(-1, 3)".
        # We can probably weed out most such reshapings.
        # Grep for the above expression in order to find places that should
        # be updated.

    @abstractmethod
    def get_value(self) -> float:
        """Return function value at current coordinates."""

    @abstractmethod
    def iterimages(self):
        """Yield domain objects that can be saved as trajectory.

        For example this can yield Atoms objects if the optimizer
        has a trajectory that can write Atoms objects."""

    def converged(self, gradient: np.ndarray, fmax: float) -> bool:
        """Standard implementation of convergence criterion.

        This assumes that forces are the actual (Nx3) forces.
        We can hopefully change this."""
        assert gradient.ndim == 1
        return self.gradient_norm(gradient) < fmax

    def gradient_norm(self, gradient):
        forces = gradient.reshape(-1, 3)  # XXX Cartesian
        return np.linalg.norm(forces, axis=1).max()

    def __ase_optimizable__(self) -> 'Optimizable':
        """Return self, being already an Optimizable."""
        return self
