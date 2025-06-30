import numpy as np
import pytest

from ase.utils.abc import Optimizable


class BoothFunctionOptimizable(Optimizable):
    """

    https://en.wikipedia.org/wiki/Test_functions_for_optimization
    """

    def __init__(self, x0):
        self.xy = np.array(x0)

    def get_x(self):
        return np.array([*self.xy, 0])  # XXX

    def set_x(self, x):
        self.xy[:] = x[:2]  # XXX

    @property
    def x(self):
        return self.xy[0]

    @property
    def y(self):
        return self.xy[1]

    @staticmethod
    def ab(x, y):
        return x + 2 * y - 7, 2 * x + y - 5

    def get_value(self):
        a, b = self.ab(*self.xy)
        return a * a + b * b

    def get_gradient(self):
        x, y = self.xy
        a, b = self.ab(*self.xy)
        # XXX negative gradient
        return -np.array([2 * a + 4 * b, 4 * a + 2 * b, 0])  # XXX

    def iterimages(self):
        return iter([])

    def ndofs(self):
        return 3  # XXX len(self.xy) + 1


def test_booth():
    from ase.optimize.bfgs import BFGS

    x0 = [1.234, 2.345]
    target = BoothFunctionOptimizable(x0)
    opt = BFGS(target)
    opt.run(fmax=1e-9)

    assert target.xy == pytest.approx([1, 3], abs=1e-7)
    assert target.get_value() == pytest.approx(0)
