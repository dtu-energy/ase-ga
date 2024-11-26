import numpy as np
import pytest

from ase import Atom, Atoms
from ase.io.nwchem import write_nwchem_in


@pytest.fixture()
def atomic_configuration():
    molecule = Atoms(pbc=False)
    molecule.append(Atom('C', [0, 0, 0]))
    molecule.append(Atom('O', [1.6, 0, 0]))
    return molecule


@pytest.fixture()
def calculator_parameters():
    params = dict(memory='1024 mb',
                  dft=dict(xc='b3lyp',
                           mult=1,
                           maxiter=300),
                  basis='6-311G*')
    return params


def test_echo(atomic_configuration, calculator_parameters, tmpdir):
    fd = tmpdir.mkdir('sub').join('nwchem.in')
    write_nwchem_in(
        fd,
        atomic_configuration,
        echo=False,
        **calculator_parameters)
    content = [line.rstrip('\n') for line in fd.readlines()]
    assert 'echo' not in content

    write_nwchem_in(
        fd,
        atomic_configuration,
        echo=True,
        **calculator_parameters)
    content = [line.rstrip('\n') for line in fd.readlines()]
    assert 'echo' in content


def test_params(atomic_configuration, calculator_parameters, tmpdir):
    fd = tmpdir.mkdir('sub').join('nwchem.in')
    write_nwchem_in(fd, atomic_configuration, **calculator_parameters)
    content = [line.rstrip('\n') for line in fd.readlines()]

    for key, value in calculator_parameters.items():
        for line in content:
            flds = line.split()
            if len(flds) == 0:
                continue
            if flds[0] == key:
                break
        else:
            assert False
        if key == 'basis':  # special case
            pass
        elif isinstance(value, str):
            assert len(value.split()) == len(flds[1:])
            assert all(v == f for v, f in zip(value.split(), flds[1:]))
        elif isinstance(value, (int, float)):
            assert len(flds[1:]) == 1
            assert np.isclose(value, float(flds[1]))
