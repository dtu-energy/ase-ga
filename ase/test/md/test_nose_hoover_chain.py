from __future__ import annotations

import numpy as np
import pytest

import ase.build
import ase.units
from ase.md.nose_hoover_chain import (
    NoseHooverChainNVT,
    NoseHooverChainThermostat,
)
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary


@pytest.mark.parametrize("tchain", [1, 3])
def test_thermostat(tchain: int):
    atoms = ase.build.bulk(
        "Cu", crystalstructure='hcp', a=2.53, c=4.11
    ).repeat(2)

    timestep = 1.0 * ase.units.fs
    thermostat = NoseHooverChainThermostat(
        num_atoms_global=len(atoms),
        masses=atoms.get_masses()[:, None],
        temperature_K=1000,
        tdamp=100 * timestep,
        tchain=tchain,
    )

    rng = np.random.default_rng(0)
    p = rng.standard_normal(size=(len(atoms), 3))

    n = 1000
    p_start = p.copy()
    eta_start = thermostat._eta.copy()
    p_eta_start = thermostat._p_eta.copy()
    for _ in range(n):
        p = thermostat.integrate_nhc(p, timestep)
    for _ in range(2 * n):
        p = thermostat.integrate_nhc(p, -0.5 * timestep)

    assert np.allclose(p, p_start, atol=1e-6)
    assert np.allclose(thermostat._eta, eta_start, atol=1e-6)
    assert np.allclose(thermostat._p_eta, p_eta_start, atol=1e-6)


@pytest.mark.parametrize("tchain", [1, 3])
def test_nose_hoover_chain_nvt(asap3, tchain: int):
    atoms = ase.build.bulk("Cu").repeat((2, 2, 2))
    atoms.calc = asap3.EMT()

    temperature_K = 300
    rng = np.random.default_rng(0)
    MaxwellBoltzmannDistribution(
        atoms,
        temperature_K=temperature_K, force_temp=True, rng=rng
    )
    Stationary(atoms)

    timestep = 1.0 * ase.units.fs
    md = NoseHooverChainNVT(
        atoms,
        timestep=timestep,
        temperature_K=temperature_K,
        tdamp=100 * timestep,
        tchain=tchain,
    )
    conserved_energy1 = md.get_conserved_energy()
    md.run(100)
    conserved_energy2 = md.get_conserved_energy()
    assert np.allclose(np.sum(atoms.get_momenta(), axis=0), 0.0)
    assert np.isclose(conserved_energy1, conserved_energy2, atol=1e-3)
