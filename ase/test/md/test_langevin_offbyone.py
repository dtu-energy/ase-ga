import numpy as np
import pytest

from ase.build import bulk
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.units import fs, kB
import time

def make_atoms(T, rng, asap3):
    atoms = bulk('Cu', cubic=True)
    MaxwellBoltzmannDistribution(atoms, temperature_K=T, rng=rng)
    Stationary(atoms)
    atoms.calc = asap3.EMT()
    return atoms


class EnergyMeasurer:
    def __init__(self, atoms):
        self.atoms = atoms
        self.energies = []
        self.temperatures = []

    def __call__(self):
        energy = self.atoms.get_potential_energy() \
            + self.atoms.get_kinetic_energy()
        temperature = self.atoms.get_temperature()
        self.energies.append(energy)
        self.temperatures.append(temperature)


def run_nvt(atoms, nsteps, dt, expected_T, dynmaker, rng=None, intval=5):
    '''Run NVT dynamics, testing the behaviour.

    Parameters:
    atoms:      The atoms object.
    nsteps:     Length of simulation.
    dt:         Time step.
    expected_T: Expected temperature
    dynmaker:   Function making the dynamics object.
    rng:        Random number generator.
    intval:     Interval for taking data.
    '''
    runtime = nsteps * dt
    tau = runtime / 5    # Energy relaxation time in ideal gas

    # We pass *half* the energy relaxation time to the generator of the
    # dynamics, as we have a solid, where the relaxation time will be twice
    # that of the ideal gas, since the same amount of potential and kinetic
    # energy needs to be added to the system.
    dyn = dynmaker(atoms, expected_T, dt, tau / 2, rng=rng, loginterval=nsteps//10)
    measure = EnergyMeasurer(atoms)
    dyn.attach(measure, interval=1)
    dyn.run(nsteps)

    # Run again while taking data
    dyn = dynmaker(atoms, expected_T, dt, tau, rng=rng, loginterval=nsteps)
    measure = EnergyMeasurer(atoms)
    dyn.attach(measure, interval=intval)
    com_before = atoms.get_center_of_mass()
    dyn.run(nsteps * 10)
    com_after = atoms.get_center_of_mass()
    energies = measure.energies
    energies = energies[len(energies) // 10:]
    temperatures2 = measure.temperatures
    temperatures2 = temperatures2[len(temperatures2) // 5:]

    stdev_energy = np.std(energies)
    agv_temperature = np.mean(temperatures2)
    # Expected energy fluctuation: sqrt(k_B T^2 3 N k_B) = k_B * T * sqrt(3 * N)
    expected = kB * expected_T * np.sqrt(3 * len(atoms))

    # Output results
    print(f'Observed energy fluctuation: {stdev_energy:.2f} eV')
    print(f'Expected energy fluctuation: {expected:.2f} eV')
    print(f'Error: {(stdev_energy / expected - 1) * 100:.1f}%')
    #assert np.abs(stdev_energy - expected) < 0.25 * expected, \
    #    'Energy fluctuations'

    # Temperature error: We should be able to detect a error of 1/N_atoms
    # The factor .67 is arbitrary, smaller than 1.0 so we consistently
    # detect errors, but not so small that we get false positives.
    maxtemperr = 0.67 * 1 / len(atoms)
    # ... but not if we don't have good statistics.
    print(f'Observed average temperature:  {agv_temperature:.2f} K'
              + f'   (expected {expected_T:.2f} K)')
    print(f'Error: {(agv_temperature / expected_T - 1) * 100:.1f}%  '
              + f'(max: {maxtemperr * 100:.1f}%)')
    assert np.abs(agv_temperature - expected_T) < expected_T * maxtemperr, \
        'Temperature'

    print('Center of mass before:', com_before)
    print('Center of mass after: ', com_after)


def make_langevin(atoms, desired_T, dt, tau, rng, loginterval):
    # tau is the energy relaxation time.  The velocity relaxation time
    # should be the double.
    return Langevin(atoms,
                    dt,
                    temperature_K=desired_T,
                    friction=1 / (2 * tau),
                    logfile='-',
                    loginterval=loginterval,
                    rng=rng,
    )


@pytest.mark.slow()
def test_langevin_offbyone(asap3):
    T0 = 300
    dt = 5 * fs
    nsteps = 1000
    rng = np.random.default_rng(2718281828459045)

    runtime = time.perf_counter()

    atoms = make_atoms(T0 / 10, rng, asap3)
    run_nvt(atoms, nsteps, dt, T0, dynmaker=make_langevin,
            rng=rng)

    runtime = time.perf_counter() - runtime
    print(f'Runtime: {runtime:.2f} s.')