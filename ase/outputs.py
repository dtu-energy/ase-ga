"""Module for ``Property`` and ``Properties``."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping
from typing import Sequence, Union

import numpy as np


class Properties(Mapping):
    def __init__(self, dct: dict) -> None:
        self._dct: dict[str, Property] = {}
        for name, value in dct.items():
            self._setvalue(name, value)

    def __len__(self) -> int:
        return len(self._dct)

    def __iter__(self):
        return iter(self._dct)

    def __getitem__(self, name) -> Property:
        return self._dct[name]

    def _setvalue(self, name: str, value) -> None:
        if name in self._dct:
            # Which error should we raise for already existing property?
            raise ValueError(f'{name} already set')

        prop = all_outputs[name]
        value = prop.normalize_type(value)
        shape = np.shape(value)

        if not self.shape_is_consistent(prop, value):
            raise ValueError(f'{name} has bad shape: {shape}')

        for i, spec in enumerate(prop.shapespec):
            if not isinstance(spec, str) or spec in self._dct:
                continue
            self._setvalue(spec, shape[i])

        self._dct[name] = value

    def shape_is_consistent(self, prop: Property, value) -> bool:
        """Return whether shape of values is consistent with properties.

        For example, forces of shape (7, 3) are consistent
        unless properties already have "natoms" with non-7 value.
        """
        shapespec = prop.shapespec
        shape = np.shape(value)
        if len(shapespec) != len(shape):
            return False
        for dimspec, dim in zip(shapespec, shape):
            if isinstance(dimspec, str):
                dimspec = self._dct.get(dimspec, dim)
            if dimspec != dim:
                return False
        return True

    def __repr__(self) -> str:
        clsname = type(self).__name__
        return f'({clsname}({self._dct})'


all_outputs: dict[str, Property] = {}


class Property(ABC):
    def __init__(self, name: str, dtype: type, shapespec: tuple) -> None:
        self.name = name
        if dtype not in {float, int}:  # Others?
            raise ValueError(dtype)
        self.dtype = dtype
        self.shapespec = shapespec

    @abstractmethod
    def normalize_type(self, value): ...

    def __repr__(self) -> str:
        typename = self.dtype.__name__  # Extend to other than float/int?
        shape = ', '.join(str(dim) for dim in self.shapespec)
        return f'Property({self.name!r}, dtype={typename}, shape=[{shape}])'


class ScalarProperty(Property):
    def __init__(self, name: str, dtype: type) -> None:
        super().__init__(name, dtype, ())

    def normalize_type(self, value):
        if not np.isscalar(value):
            raise TypeError('Expected scalar')
        return self.dtype(value)


class ArrayProperty(Property):
    def normalize_type(self, value):
        if np.isscalar(value):
            raise TypeError('Expected array, got scalar')
        return np.asarray(value, dtype=self.dtype)


ShapeSpec = Union[str, int]


def _defineprop(
    name: str,
    dtype: type = float,
    shape: Union[ShapeSpec, Sequence[ShapeSpec]] = (),
) -> Property:
    """Create, register, and return a property."""

    if isinstance(shape, (int, str)):
        shape = (shape,)

    shape = tuple(shape)
    prop: Property
    if len(shape) == 0:
        prop = ScalarProperty(name, dtype)
    else:
        prop = ArrayProperty(name, dtype, shape)

    assert name not in all_outputs, name
    all_outputs[name] = prop
    return prop


# Atoms, energy, forces, stress:
_defineprop('natoms', int)
_defineprop('energy', float)
_defineprop('energies', float, shape='natoms')
_defineprop('free_energy', float)
_defineprop('forces', float, shape=('natoms', 3))
_defineprop('stress', float, shape=6)
_defineprop('stresses', float, shape=('natoms', 6))

# Electronic structure:
_defineprop('nbands', int)
_defineprop('nkpts', int)
_defineprop('nspins', int)
_defineprop('fermi_level', float)
_defineprop('kpoint_weights', float, shape='nkpts')
_defineprop('ibz_kpoints', float, shape=('nkpts', 3))
_defineprop('eigenvalues', float, shape=('nspins', 'nkpts', 'nbands'))
_defineprop('occupations', float, shape=('nspins', 'nkpts', 'nbands'))

# Miscellaneous:
_defineprop('dipole', float, shape=3)
_defineprop('charges', float, shape='natoms')
_defineprop('magmom', float)
_defineprop('magmoms', float, shape='natoms')  # XXX spinors?
_defineprop('polarization', float, shape=3)
_defineprop('dielectric_tensor', float, shape=(3, 3))
_defineprop('born_effective_charges', float, shape=('natoms', 3, 3))

# We might want to allow properties that are part of Atoms, such as
# positions, numbers, pbc, cell.  It would be reasonable for those
# concepts to have a formalization outside the Atoms class.


# def to_singlepoint(self, atoms):
#    from ase.calculators.singlepoint import SinglePointDFTCalculator
#    return SinglePointDFTCalculator(atoms,
#                                    efermi=self.fermi_level,

# We can also retrieve (P)DOS and band structure.  However:
#
# * Band structure may require bandpath, which is an input, and
#   may not necessarily be easy or possible to reconstruct from
#   the outputs.
#
# * Some calculators may produce the whole BandStructure object in
#   one go (e.g. while parsing)
#
# * What about HOMO/LUMO?  Can be obtained from
#   eigenvalues/occupations, but some codes provide real data.  We
#   probably need to distinguish between HOMO/LUMO inferred by us
#   versus values provided within the output.
#
# * HOMO is sometimes used as alternative reference energy for
#   band structure.
#
# * What about spin-dependent (double) Fermi level?
#
# * What about 3D arrays?  We will almost certainly want to be
#   connected to an object that can load dynamically from a file.
