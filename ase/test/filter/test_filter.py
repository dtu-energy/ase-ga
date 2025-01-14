import pytest

from ase.build import molecule
from ase.calculators.emt import EMT
from ase.filters import Filter
from ase.optimize import QuasiNewton


@pytest.mark.optimize()
def test_filter(testdir):
    """Test that the filter and trajectories are playing well together."""

    atoms = molecule('CO2')
    atoms.calc = EMT()
    filter = Filter(atoms, indices=[1, 2])

    with QuasiNewton(filter, trajectory='filter-test.traj',
                     logfile='filter-test.log') as opt:
        opt.run()
    # No assertions=??
