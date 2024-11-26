import numpy as np
import pytest

from ase import Atoms
from ase.build import bulk
from ase.md.npt import NPT
from ase.md.nptberendsen import NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.units import GPa, bar, fs


@pytest.fixture(scope='module')
def berendsenparams():
    """Parameters for the two Berendsen algorithms."""
    Bgold = 220.0 * GPa  # Bulk modulus of gold, in bar (1 GPa = 10000 bar)
    nvtparam = dict(temperature_K=300, taut=1000 * fs)
    nptparam = dict(temperature_K=300, pressure_au=5000 * bar, taut=1000 * fs,
                    taup=1000 * fs,
                    compressibility_au=1 / Bgold)
    return dict(nvt=nvtparam, npt=nptparam)


@pytest.fixture(scope='module')
def equilibrated(asap3, berendsenparams):
    """Make an atomic system with equilibrated temperature and pressure."""
    rng = np.random.RandomState(42)
    # Must be big enough to avoid ridiculous fluctuations
    atoms = bulk('Au', cubic=True).repeat((3, 3, 3))
    atoms.calc = asap3.EMT()
    MaxwellBoltzmannDistribution(atoms, temperature_K=100, force_temp=True,
                                 rng=rng)
    Stationary(atoms)
    assert abs(atoms.get_temperature() - 100) < 0.0001
    with NPTBerendsen(atoms, timestep=20 * fs, logfile='-',
                      loginterval=200,
                      **berendsenparams['npt']) as md:
        # Equilibrate for 20 ps
        md.run(steps=1000)
    T = atoms.get_temperature()
    pres = -atoms.get_stress(
        include_ideal_gas=True)[:3].sum() / 3 / GPa * 10000
    print(f"Temperature: {T:.2f} K    Pressure: {pres:.2f} bar")
    return atoms


def propagate(atoms, asap3, algorithm, algoargs):
    T = []
    p = []
    with algorithm(
            atoms,
            timestep=20 * fs,
            logfile='-',
            loginterval=1000,
            **algoargs) as md:
        # Gather data for 50 ps
        for _ in range(500):
            md.run(5)
            T.append(atoms.get_temperature())
            pres = - atoms.get_stress(include_ideal_gas=True)[:3].sum() / 3
            p.append(pres)
    Tmean = np.mean(T)
    p = np.array(p)
    pmean = np.mean(p)
    print('Temperature: {:.2f} K +/- {:.2f} K  (N={})'.format(
        Tmean, np.std(T), len(T)))
    print('Center-of-mass corrected temperature: {:.2f} K'.format(
        Tmean * len(atoms) / (len(atoms) - 1)))
    print('Pressure: {:.2f} bar +/- {:.2f} bar  (N={})'.format(
        pmean / bar, np.std(p) / bar, len(p)))
    return Tmean, pmean


# Not a real optimizer test but uses optimizers.
# We should probably not mark this (in general)
@pytest.mark.optimize()
@pytest.mark.slow()
def test_nvtberendsen(asap3, equilibrated, berendsenparams, allraise):
    t, _ = propagate(Atoms(equilibrated), asap3,
                     NVTBerendsen, berendsenparams['nvt'])
    assert abs(t - berendsenparams['nvt']['temperature_K']) < 0.5


@pytest.mark.optimize()
@pytest.mark.slow()
def test_nptberendsen(asap3, equilibrated, berendsenparams, allraise):
    t, p = propagate(Atoms(equilibrated), asap3,
                     NPTBerendsen, berendsenparams['npt'])
    assert abs(t - berendsenparams['npt']['temperature_K']) < 1.0
    assert abs(p - berendsenparams['npt']['pressure_au']) < 25.0 * bar


@pytest.mark.optimize()
@pytest.mark.slow()
def test_npt(asap3, equilibrated, berendsenparams, allraise):
    params = berendsenparams['npt']
    # NPT uses different units.  The factor 1.3 is the bulk modulus of gold in
    # ev/Å^3
    t, p = propagate(Atoms(equilibrated), asap3, NPT,
                     dict(temperature_K=params['temperature_K'],
                          externalstress=params['pressure_au'],
                          ttime=params['taut'],
                          pfactor=params['taup']**2 * 1.3))
    # Unlike NPTBerendsen, NPT assumes that the center of mass is not
    # thermalized, so the kinetic energy should be 3/2 ' kB * (N-1) * T
    n = len(equilibrated)
    assert abs(t - (n - 1) / n * berendsenparams['npt']['temperature_K']) < 1.0
    assert abs(p - berendsenparams['npt']['pressure_au']) < 100.0 * bar
