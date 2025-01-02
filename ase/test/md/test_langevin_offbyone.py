import numpy as np
import pytest

# matplotlib only imported in debugging mode (in function test_nvt)
# import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from ase.build import bulk
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.units import fs, kB


def make_atoms(T, rng, asap3):
    atoms = bulk('Cu', cubic=True)
    MaxwellBoltzmannDistribution(atoms, temperature_K=T, rng=rng)
    Stationary(atoms)
    atoms.calc = asap3.EMT()
    return atoms


class MeasureEnergy:
    def __init__(self, atoms):
        self.atoms = atoms
        self.energies = []
        self.temperatures = []

    def __call__(self):
        e = self.atoms.get_potential_energy() + self.atoms.get_kinetic_energy()
        T = self.atoms.get_temperature()
        self.energies.append(e)
        self.temperatures.append(T)


def tempcurve(t, A, tau, B):
    'Temperature curve to fit to'
    return A * np.exp(-t / tau) + B


def run_nvt(atoms, nsteps, dt, T0, dynmaker, rng=None,
            intval=5, plot=False, sloppytime=False, failfluct=False):
    '''Run NVT dynamics, testing the behaviour.

    Parameters:
    atoms:      The atoms object.
    nsteps:     Length of simulation.
    dt:         Time step.
    T0:         Expected temperature
    dynmaker:   Function making the dynamics object.
    rng:        Random number generator.
    intval:     Interval for taking data.
    plot:       Plot the swing-in temperature graph (default: False).
    sloppytime: Test of swing-in time is sloppy (+/- 50% instead of +/- 1%).
                (for nondeterministic dynamics)
    failfluct:  Test of energy fluctuations is expected to fail.
                (for dynamics that does not produce a true Canonical Ensemble)
    '''
    runtime = nsteps * dt
    tau = runtime / 10    # Energy relaxation time in ideal gas

    # We pass *half* the energy relaxation time to the generator of the
    # dynamics, as we have a solid, where the relaxation time will be twice
    # that of the ideal gas, since the same amount of potential and kinetic
    # energy needs to be added to the system.
    dyn = dynmaker(atoms, T0, dt, tau / 2, rng=rng, logint=250)
    measure = MeasureEnergy(atoms)
    dyn.attach(measure, interval=1)
    dyn.run(nsteps)
    temperatures = measure.temperatures
    del dyn

    # Fit temperature curve
    temperatures = temperatures[5:]
    times = np.arange(len(temperatures)) * dt
    (DeltaTfit, tau_fit, T_fit), _ = curve_fit(tempcurve, times, temperatures,
                                                   (-T0, tau, T0))

    # Run again with smaller tau
    tausmall = runtime / 10
    dyn = dynmaker(atoms, T0, dt, tausmall, rng=rng, logint=2500)
    measure = MeasureEnergy(atoms)
    dyn.attach(measure, interval=intval)
    com_before = atoms.get_center_of_mass()
    if failfluct:
        # No need for good statistics if it fails anyway
        dyn.run(nsteps)
    else:
        dyn.run(nsteps * 10)
    com_after = atoms.get_center_of_mass()
    energies = measure.energies
    energies = energies[len(energies) // 10:]
    temperatures2 = measure.temperatures
    temperatures2 = temperatures2[len(temperatures2) // 5:]
    times2 = np.arange(len(temperatures2)) * dt * intval
    del dyn

    stdE = np.std(energies)
    avgT = np.mean(temperatures2)
    # Expected energy fluctuation: sqrt(k_B T^2 3 N k_B) = k_B * T * sqrt(3 * N)
    expected = kB * T0 * np.sqrt(3 * len(atoms))

    # Output results
    print(f'Part 1 temperature:   {T_fit:.2f} K    (expected {T0:.2f})')
    print(f'Time constant:       {tau_fit / fs:.1f} fs  '
              + '(expected {tau / fs:.1f}  '
              + 'error {(tau_fit / tau - 1) * 100:.1f}%)')
    print(f'Initial temperature offset:  {DeltaTfit:.2f} K')
    print('Note: Due to the small system size, the numbers above will')
    print('      be completely off.')
    print()
    print(f'Observed energy fluctuation: {stdE:.2f} eV')
    print(f'Expected energy fluctuation: {expected:.2f} eV')
    print(f'Error: {(stdE / expected - 1) * 100:.1f}%')
    if failfluct:
        print('Failed fluctuations EXPECTED for this dynamics!')
    else:
        assert np.abs(stdE - expected) < 0.25 * expected, 'Energy fluctuations'

    # Temperature error: We should be able to detect a error of 1/N_atoms
    maxtemperr = 2 / 3 * 3 / atoms.get_number_of_degrees_of_freedom()
    # ... but not if we don't have good statistics.
    if failfluct:
        maxtemperr *= 3
    print(f'Observed average temperature:  {avgT:.2f} K'
              + '   (expected {T0:.2f} K)')
    print(f'Error: {(avgT / T0 - 1) * 100:.1f}%  '
              + '(max: {maxtemperr * 100:.1f}%)')
    assert np.abs(avgT - T0) < T0 * maxtemperr, 'Temperature'

    print('Center of mass before:', com_before)
    print('Center of mass after: ', com_after)

    if plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(times, temperatures, 'b.')
        ax.plot(times, tempcurve(times, DeltaTfit, tau_fit, T_fit), 'k-')
        fig2, ax2 = plt.subplots()
        ax2.plot(times2, temperatures2, 'b.')
        plt.show(block=True)
        del fig, fig2


def lgvdynmaker(atoms, T0, dt, tau, rng, logint):
    # tau is the energy relaxation time.  The velocity relaxation time
    # should be the double.
    return Langevin(atoms, dt, temperature_K=T0, friction=1 / (2 * tau),
                        logfile='-', loginterval=logint, rng=rng)


@pytest.mark.slow()
def test_langevin_offbyone(asap3):
    T0 = 300
    dt = 5 * fs
    nsteps = 10000
    rng = np.random.default_rng(2718281828459045)

    atoms = make_atoms(T0 / 10, rng, asap3)
    run_nvt(atoms, nsteps, dt, T0, dynmaker=lgvdynmaker,
            rng=rng, plot=False, sloppytime=True)
