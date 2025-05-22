.. module:: ase.calculators.harmonic

.. _harmonic:

===================
Harmonic calculator
===================

Introduction
============

The local Harmonic Approximation of the potential energy surface (PES) is
commonly applied in atomistic simulations to estimate entropy, i.e. free
energy, at elevated temperatures (e.g. in ASE via :mod:`~ase.thermochemistry`).
The term 'harmonic' refers to a second order Taylor series of the PES for a
local reference configuration in Cartesian coordinates expressed in a Hessian
matrix. With the Hessian matrix (e.g. computed numerically in ASE via
:mod:`~ase.vibrations`) normal modes and harmonic vibrational frequencies can
be obtained.

The following :class:`HarmonicCalculator` can be used to compute energy and forces
with a Hessian-based harmonic force field (:class:`HarmonicForceField`).
Moreover, it can be used to compute Anharmonic Corrections to the
Harmonic Approximation. [1]_

.. [1] Amsler, J. et al., Anharmonic Correction to Adsorption Free Energy
       from DFT-Based MD Using Thermodynamic Integration,
       J. Chem. Theory Comput. 2021, 17 (2), 1155-1169.
       :doi:`10.1021/acs.jctc.0c01022`.

.. autoclass:: ase.calculators.harmonic.HarmonicCalculator

.. autoclass:: ase.calculators.harmonic.HarmonicForceField

.. note::

   The reference Hessians in **x** and **q** can be inspected via
   ``HarmonicForceField.hessian_x`` and ``HarmonicForceField.hessian_q``.

Theory for Anharmonic Correction via Thermodynamic Integration (TI)
===================================================================
Thermodynamic integration (TI), i.e. `\lambda`-path integration,
connects two thermodynamic states via a `\lambda`-path.
Here, the TI begins from a reference system '0' with known free energy
(Harmonic Approximation) and the Anharmonic Correction is obtained via
integration over the `\lambda`-path to the target system '1' (the fully
interacting anharmonic system).
Hence, the free energy of the target system can be written as

.. math::
    A_1 = A_0 + \Delta A_{0 \rightarrow 1}

where the second term corresponds to the integral over the `\lambda`-path

.. math::

    \Delta A_{0 \rightarrow 1} = \int_0^1 d \lambda
    \langle H_1 - H_0 \rangle_\lambda

The term `\langle ... \rangle_\lambda` represents the NVT ensemble
average of the system driven by the classical Hamiltonian
`\mathcal{H}_\lambda` determined by the coupling parameter
`\lambda \in [0,1]`

.. math::

    \mathcal{H}_\lambda = \lambda \mathcal{H}_1 + (1 - \lambda) \mathcal{H}_0

Since the Hamiltonians differ only in their potential energy contributions
`V_1` and `V_0`, the free energy change can be computed from the
potentials

.. math::

    \Delta A_{0 \rightarrow 1} = \int_0^1 d \lambda
    \langle V_1 - V_0 \rangle_\lambda

The Cartesian coordinates **x** used in the common Harmonic Approximation are
not insensitive to overall rotations and translations that must leave the total
energy invariant.
This limitation can be overcome by transformation of the Hessian in **x**
to a suitable coordinate system **q** (e.g. internal coordinates).
Since the force field of that Hessian which is harmonic in **x** is not
necessarily equivalently harmonic in **q**, the free energy correction can be
rewritten to

.. math::
    A_1 = A_{0,\mathbf{x}} + \Delta A_{0,\mathbf{x} \rightarrow 0,\mathbf{q}}
    + \Delta A_{0,\mathbf{q} \rightarrow 1}

The terms in this equation correspond to the free energy from the Harmonic
Approximation with the reference Hessian (`A_{0,\mathbf{x}}`), the free
energy change due to the coordinate transformation
(`\Delta A_{0,\mathbf{x} \rightarrow 0,\mathbf{q}}`) obtained via TI
(see Example 3) and the free energy change from the harmonic to the fully
interacting system (`\Delta A_{0,\mathbf{q} \rightarrow 1}`) obtained via
TI (see Example 4).
Please see Amsler, J. et al. for details. [1]_

.. note::

    Anharmonicity is quantified by comparison of the total free energy
    `A_1` to the free energy contributions by the standard Harmonic
    Approximation with the unmodified Hessian.
    The reference Hessian and its free energy contribution
    `A_{0,\mathbf{x}}` have no meaning outside the TI procedure.

Examples
========
Prerequisites: :class:`~ase.Atoms` object (``ref_atoms``),
its energy (``ref_energy``) and Hessian (``hessian_x``).

Example 1: Cartesian coordinatates
----------------------------------
In Cartesian coordinates, forces and energy are not invariant with respect
to rotations and translations of the system.

.. code-block:: python

    import numpy as np
    from ase.calculators.harmonic import HarmonicForceField, HarmonicCalculator
    hff = HarmonicForceField(ref_atoms=ref_atoms, ref_energy=ref_energy,
                             hessian_x=hessian_x)
    atoms = ref_atoms.copy()
    atoms.calc = HarmonicCalculator(hff)

.. note::

   Forces and energy can be computed via :meth:`~ase.Atoms.get_forces` and
   :meth:`~ase.Atoms.get_potential_energy` for any configuration that does
   not involve rotations with respect to the configuration of ``ref_atoms``.

.. warning::

   In case of system rotations, Cartesian coordinates return incorrect values
   and thus cannot be used without an appropriate coordinate system
   as demonstrated in the Supporting Information of Amsler, J. et al.. [1]_

Example 2: Internal Coordinates
-------------------------------
To compute forces and energy correctly even for rotated systems,
a user-defined coordinate system must be provided.
Within this coordinate system, energy and forces must be invariant with
respect to rotations and translations of the system.
For this purpose internal coordinates (distances, angles, dihedrals,
coordination numbers and linear combinations thereof, etc.) are widely used.
The following example works on a water molecule (:mol:`H_2O`) stored in
``ref_atoms``.

.. literalinclude:: ../../../ase/test/calculator/test_harmonic.py
    :language: python
    :start-after: start doc example 3
    :end-before: end doc example 3

.. literalinclude:: ../../../ase/test/calculator/test_harmonic.py
    :language: python
    :start-after: test_internals():
    :end-before: atoms = setup_water(calc)  # distorted copy of ref_atoms
    :dedent: 4

Example 3: Free Energy Change due to Coordinate Transformation
--------------------------------------------------------------
A transformation of the coordinate system may transform the force field.
The change in free energy due to this transformation
(`\Delta A_{0,\mathbf{x} \rightarrow 0,\mathbf{q}}`) can be computed via
thermodynamic (`\lambda`-path) integration. [1]_

.. literalinclude:: ../../../ase/test/calculator/test_harmonic.py
    :language: python
    :start-after: test_thermodynamic_integration():
    :end-before: assert -0.005 < dA < 0.005
    :dedent: 4

Integration of the mean energy differences ('ediffs') over the integration grid
(`\lambda` path) leads to the change in free energy due to the coordinate
transformation.

Example 4: Anharmonic Corrections
---------------------------------
The strategy in Example 3 can be used to compute anharmonic corrections to the
Harmonic Approximation when the :class:`HarmonicCalculator` is coupled with
a calculator that can compute interactions beyond the Harmonic Approximation,
e.g. :mod:`~ase.calculators.vasp`.

.. note::

   The obtained Anharmonic Correction applies to the Harmonic Approximation
   (`A_{0,\mathbf{x}}`) of the reference system with the reference Hessian which
   is generated during initialization of the Calculator and
   may differ from the standard Harmonic Approximation.
   The vibrations for the reference system can be computed numerically with
   high accuracy.

    >>> from ase.vibrations import Vibrations
    >>> atoms = ref_atoms.copy()
    >>> atoms.calc = calc_harmonic_0  # with cartesian=True
    >>> vib = Vibrations(atoms, nfree=4, delta=1e-5)
    >>> vib.run()
