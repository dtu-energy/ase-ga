"""Tools for locally structure optimization."""

import os
from time import time

import numpy as np

from ase.build import niggli_reduce
from ase.calculators.calculator import all_changes
from ase.calculators.lj import LennardJones
from ase.calculators.singlepoint import SinglePointCalculator
from ase.ga import set_raw_score
from ase.optimize.precon import PreconLBFGS
from ase.units import kB


def relax(atoms):
    """Performs a variable-cell optimization of the Atoms object."""
    t1 = time()

    # LJ parameters from G. Galassi and D.J. Tildesley,
    # "Phase Diagrams of Diatomic Molecules:
    #  Using the Gibbs Ensemble Monte Carlo Method",
    # Molecular Simulations, 13, 11 (1994).
    calc = HarmonicPlusLennardJones(
        epsilon=37.3 * kB, sigma=3.31, rc=12.0, r0=1.12998, k=10.0
    )
    atoms.calc = calc

    dyn = PreconLBFGS(
        atoms,
        variable_cell=True,
        maxstep=0.2,
        use_armijo=True,
        logfile='opt.log',
        trajectory='opt.traj',
    )
    dyn.run(fmax=3e-2, smax=5e-4, steps=250)

    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    s = atoms.get_stress()
    finalize(atoms, energy=e, forces=f, stress=s)

    os.system('mv opt.log prev.log')
    os.system('mv opt.traj prev.traj')

    t2 = time()
    print('Relaxation took %.3f seconds.' % (t2 - t1), flush=True)


def finalize(atoms, energy=None, forces=None, stress=None):
    niggli_reduce(atoms)
    atoms.wrap()
    calc = SinglePointCalculator(
        atoms, energy=energy, forces=forces, stress=stress
    )
    atoms.calc = calc
    raw_score = -atoms.get_potential_energy()
    set_raw_score(atoms, raw_score)


class HarmonicPlusLennardJones(LennardJones):
    """Lennard-Jones potential for intermolecular interactions
    (parameters: epsilon, sigma, rc) and a harmonic potential
    for intramolecular interactions (parameters: k, r0).

    Only works for structures consisting of a series
    of molecular dimers and with only one element.
    """

    implemented_properties = ['energy', 'forces', 'stress']
    default_parameters = {'k': 1.0, 'r0': 1.0}
    nolabel = True

    def __init__(self, **kwargs):
        LennardJones.__init__(self, **kwargs)

    def calculate(
        self, atoms=None, properties=['energy'], system_changes=all_changes
    ):
        LennardJones.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)

        epsilon = self.parameters.epsilon
        sigma = self.parameters.sigma
        k = self.parameters.k
        r0 = self.parameters.r0
        rc = self.parameters.rc
        if rc is None:
            rc = 3 * r0

        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((3, 3))

        tags = self.atoms.get_tags()
        for tag in np.unique(tags):
            # Adding (intramolecular) harmonic potential part
            indices = np.where(tags == tag)[0]
            assert len(indices) == 2
            a1, a2 = indices
            d = self.atoms.get_distance(a1, a2, mic=True, vector=True)
            r = np.linalg.norm(d)
            energy += 0.5 * k * (r - r0) ** 2
            f = -k * (r - r0) * d
            forces[a1] -= f
            forces[a2] += f
            stress += np.dot(np.array([f]).T, np.array([d]))

            # Substracting intramolecular LJ part
            r2 = r**2
            c6 = (sigma**2 / r2) ** 3
            c12 = c6**2
            energy += -4 * epsilon * (c12 - c6).sum()
            f = (24 * epsilon * (2 * c12 - c6) / r2) * d
            forces[a1] -= -f
            forces[a2] += -f
            stress += -np.dot(np.array([f]).T, np.array([d]))

        if 'stress' in properties:
            stress += stress.T.copy()
            stress *= -0.5 / self.atoms.get_volume()
            self.results['stress'] += stress.flat[[0, 4, 8, 5, 2, 1]]

        self.results['energy'] += energy
        self.results['free_energy'] += energy
        self.results['forces'] += forces
