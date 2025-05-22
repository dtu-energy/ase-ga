.. module:: ase.dft.wannier
   :synopsis: Maximally localized Wannier functions

=====================================
Maximally localized Wannier functions
=====================================

This page describes how to construct the Wannier orbitals using the
class :class:`Wannier`. The page is organized as follows:

* `Introduction`_: A short summary of the basic theory.
* `The Wannier class`_ : A description of how the Wannier class is
  used, and the methods defined within.


Introduction
============

The point of Wannier functions is the transform the extended Bloch
eigenstates of a DFT calculation, into a smaller set of states
designed to facilitate the analysis of e.g. chemical bonding. This is
achieved by designing the Wannier functions to be localized in real
space instead of energy (which would be the eigen states).

The standard Wannier transformation is a unitary rotation of the Bloch
states. This implies that the Wannier functions (WF) span the same
Hilbert space as the Bloch states, i.e. they have the same eigenvalue
spectrum, and the original Bloch states can all be exactly reproduced
from a linear combination of the WF. For maximally localized Wannier
functions (MLWF), the unitary transformation is chosen such that the
spread of the resulting WF is minimized.

The standard choice is to make a unitary transformation of the
occupied bands only, thus resulting in as many WF as there are
occupied bands. If you make a rotation using more bands, the
localization will be improved, but the number of wannier functions
increase, thus making orbital based analysis harder.

The class defined here allows for construction of *partly* occupied
MLWF. In this scheme the transformation is still a unitary rotation
for the lowest states (the *fixed space*), but it uses a dynamically
optimized linear combination of the remaining orbitals (the *active
space*) to improve localization. This implies that e.g. the
eigenvalues of the Bloch states contained in the fixed space can be
exactly reproduced by the resulting WF, whereas the largest
eigenvalues of the WF will not necessarily correspond to any "real"
eigenvalues (this is irrelevant, as the fixed space is usually chosen
large enough, i.e. high enough above the fermilevel, that the
remaining DFT eigenvalues are meaningless anyway).

The theory behind this method can be found in
"Partly Occupied Wannier Functions" Thygesen, Hansen, and Jacobsen,
*Phys. Rev. Lett*, Vol. **94**, 26405 (2005),
:doi:`10.1103/PhysRevLett.94.026405`.

An improved localization method based on penalizing the
spread functional proportionally to the variance of spread distributions
has since been published as
"Spread-balanced Wannier functions:
Robust and automatable orbital localization"
Fontana, Larsen, Olsen, and Thygesen,
*Phys. Rev. B* **104**, 125140 (2021),
:doi:`10.1103/PhysRevB.104.125140`.


The Wannier class
=================

Usual invocation::

  from ase.dft import Wannier
  wan = Wannier(nwannier=18, calc=GPAW('save.gpw'), fixedstates=15)
  wan.localize() # Optimize rotation to give maximal localization
  wan.save('wannier.json') # Save localization and rotation matrix

  # Re-load using saved wannier data
  wan = Wannier(nwannier=18, calc=calc, fixedstates=15, file='wannier.json')

  # Write a cube file
  wan.write_cube(index=5, fname='wannierfunction5.cube')

For examples of how to use the **Wannier** class, see the
:ref:`wannier tutorial` tutorial.

To enable the improved localization method from the 2021 paper,
use ``Wannier(functional='var', ...)``.
Computations from the 2021 paper were performed using
``Wannier(functional='var', initialwannier='orbitals', ...)``.

.. autoclass:: Wannier
   :members:

.. note:: For calculations using **k**-points, make sure that the
   `\Gamma`-point is included in the **k**-point grid.
   The Wannier module does not support **k**-point reduction by symmetry, so
   you must use the ``usesymm=False`` keyword in the calc, and
   shift all **k**-points by a small amount (but not less than 2e-5
   in) in e.g. the x direction, before performing the calculation.
   If this is not done the symmetry program will still use time-reversal
   symmetry to reduce the number of **k**-points by a factor 2.
   The shift can be performed like this::

     from ase.dft.kpoints import monkhorst_pack
     kpts = monkhorst_pack((15, 9, 9)) + [2e-5, 0, 0]
