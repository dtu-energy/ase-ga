# fmt: off

import collections
from abc import abstractmethod

import numpy as np

# Due to the high prevalence of cyclic imports surrounding ase.optimize,
# we define the Optimizable ABC here in utils.
# Can we find a better way?


class Optimizable(collections.abc.Sized):
    @abstractmethod
    def get_positions(self):
        ...

    def get_x(self):
        return self.get_positions().ravel()

    @abstractmethod
    def set_positions(self, positions):
        ...

    def set_x(self, x):
        self.set_positions(x.reshape(-1, 3))

    @abstractmethod
    def get_gradient(self):
        # Callers who want Nx3 will do ".get_gradient().reshape(-1, 3)".
        # We can probably weed out most such reshapings.
        # Grep for the above expression in order to find places that should
        # be updated.
        ...

    @abstractmethod
    def get_value(self):
        ...

    @abstractmethod
    def iterimages(self):
        ...

    def converged(self, forces, fmax):
        return np.linalg.norm(forces, axis=1).max() < fmax

    def __ase_optimizable__(self):
        return self
