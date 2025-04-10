.. A new scriv changelog fragment.
..

I/O
---

 - :mod:`ase.db`: Unique IDs are now based on UUID rather than pseudorandom numbers that could become equal due to seeding (:mr:`3614`).
 - :mod:`ase.db`: Fix bug where unique_id could be converted to float or int (:mr:`3613`).
 - Vasp: More robust reading of CHGCAR (:mr:`3607`).
 - Lammpsdump: Read timestep from lammpsdump and set element based on mass (:mr:`3529`).
 - Vasp: Read and write velocities (:mr:`3597`).
 - DB: Support for LMDB via `ase-db-backends` project (:mr:`3564`, :mr:`3639`).
 - Espresso: Fix bug reading `alat` in some cases (:mr:`3562`).
 - GPAW: Fix reading of total charge from text file (:mr:`3519`).
 - extxyz: Somewhat restrict what properties are automatically written (:mr:`3516`).
 - Lammpsdump: Read custom property/atom LAMMPS dump data (:mr:`3510`).

Calculators
-----------

 - More robust reading of Castep XC functional (:mr:`3612`).
 - More robust saving of calculators to e.g. trajectories (:mr:`3610`).
 - Lammpslib: Fix outdated MPI check (:mr:`3594`).
 - Morse: Optionally override neighbor list implementation (:mr:`3593`).
 - EAM: Calculate stress (:mr:`3581`).

Optimizers
----------

 - Fix step counting in the
   :class:`~ase.optimize.cellawarebfgs.CellAwareBFGS` (:mr:`3588`).


Molecular dynamics
------------------

 - Improved random sampling in countour exploration (:mr:`3643`).
 - Fix small energy error in Langevin dynamics (:mr:`3567`).
 - Isotropic NPT with MTK equations (:mr:`3550`).
 - Bussi dynamics now work in parallel (:mr:`3569`).
 - Improvements to documentation (:mr:`3566`).
 - Make Nose-Hoover chain NVT faster and fix domain decomposition
   with Asap3 (:mr:`3571`).

GUI
---

 - Fix windowing bug on WSL (:mr:`3478`).


Development
-----------

 - Ruff formatter to be gradually enabled across codebase (:mr:`3600`).


Other changes
-------------

 - :meth:`~ase.cell.Cell.standard_form` can convert to upper triangular (:mr:`3623`).

 - Bugfix: :func:`~ase.geometry.geometry.get_duplicate_atoms` now respects pbc (:mr:`3609`).

 - Bugfix: Constraint masks in cell filters are now respected down to numerical precision.  Previously, the constraints could be violated by a small amount (:mr:`3603`).
 - Deprecate :func:`~ase.utils.lazyproperty` and :func:`~ase.utils.lazymethod`
   since Python now provides :func:`functools.cached_property` (:mr:`3565`).
 - Remove `nomad-upload` and `nomad-get` commands due to incompatibility
   with recent Nomad (:mr:`3563`).
 - Fix normalization of phonon DOS (:mr:`3472`).
 - :class:`~ase.io.utils.PlottingVariables` towards rotating the
   camera rather than the atoms (:mr:`2895`).

.. - A bullet item for the Other changes category.
..
.. Bugfixes
.. --------
..
.. - A bullet item for the Bugfixes category.
..
