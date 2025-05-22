.. _changelog:

=========
Changelog
=========

Git master branch
=================

.. CHANGELOG HOWTO.

   To add an entry to the changelog, create a file named
   <timestamp>_<subject>.rst inside the ase/changelog.d/ directory.
   Timestamp should be at least YYYYMMDD.

   You can also install scriv (https://pypi.org/project/scriv/) and run
   "scriv create" to do this automatically, if you do this often.

   Edit the file following a similar style to other changelog entries and
   try to choose an existing section for the release note.

   For example ase/changelog.d/20250108_amber_fix_velocities.rst with contents:

     Calculators
     -----------

     - Amber: Fix scaling of velocities in restart files (:mr:`3427`)

   For each release we generate a full changelog which is inserted below.

.. scriv-auto-changelog-start

Version 3.25.0
==============

I/O
---

- Moved Postgres, MariaDB and MySQL backends to separate project:
  https://gitlab.com/ase/ase-db-backends.  Install from PyPI with
  ``pip install ase-db-backends`` (:mr:`3545`).

- **BREAKING** ase.io.orca `read_orca_output` now returns Atoms with attached properties.
  `ase.io.read` will use this function.
  The previous behaviour (return results dictionary only) is still available from function `read_orca_outputs`. (:mr:`3599`)

- Added :func:`~ase.io.castep.write_castep_geom` and
  :func:`~ase.io.castep.write_castep_md` (:mr:`3229`)

- Fixed `:mod:`ase.data.pubchem` module to convert ``#`` in SMILES to HEX
  ``%23`` for URL (:mr:`3620`).

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

 - A new Calculator :class:`ase.calculators.tersoff.Tersoff` has been added. This is a Python implementation of a LAMMPS-style Tersoff interatomic potential. Parameters may be passed directly to the calculator as a :class:`ase.calculators.tersoff.TersoffParameters` object, or the Calculator may be constructed from a LAMMPS-style file using its ``from_lammps`` classmethod. (:mr:`3502`)

Optimizers
----------

 - Fix step counting in the
   :class:`~ase.optimize.cellawarebfgs.CellAwareBFGS` (:mr:`3588`).

 - Slightly more efficient/robust GoodOldQuasiNewton (:mr:`3570`).

Molecular dynamics
------------------

- Merged `self.communicator` into `self.comm` (:mr:`3631`).

 - Improved random sampling in countour exploration (:mr:`3643`).
 - Fix small energy error in Langevin dynamics (:mr:`3567`).
 - Isotropic NPT with MTK equations (:mr:`3550`).
 - Bussi dynamics now work in parallel (:mr:`3569`).
 - Improvements to documentation (:mr:`3566`).
 - Make Nose-Hoover chain NVT faster and fix domain decomposition
   with Asap3 (:mr:`3571`).

 - NPT now works with cells that are upper or lower triangular matrices
   (:mr:`3277`) aside from upper-only as before.

 - Fix inconsistent :meth:`irun` for NPT (:mr:`3598`).

GUI
---

 - Fix windowing bug on WSL (:mr:`3478`).

 - Added button to wrap atoms into cell (:mr:`3587`).

Development
-----------

- Changelog is now generated using ``scriv`` (:mr:`3572`).

- CI cleanup; pypi dependencies in CI jobs are now cached
  (:mr:`3628`, :mr:`3629`).
- Maximum automatic pytest workers reduced to 8 (:mr:`3628`).

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

.. scriv-auto-changelog-end


Version 3.24.0
==============

Requirements
------------

* The minimum supported Python version has increased to 3.9 (:mr:`3473`)
* Support numpy 2 (:mr:`3398`, :mr:`3400`, :mr:`3402`)
* Support spglib 2.5.0 (:mr:`3452`)

Atoms
-----
* New method :func:`~ase.Atoms.get_number_of_degrees_of_freedom()` (:mr:`3380`)
* New methods :func:`~ase.Atoms.get_kinetic_stress()`, :func:`~ase.Atoms.get_kinetic_stresses()` (:mr:`3362`)
* Prevent truncation when printing Atoms objects with 1000 or more atoms (:mr:`2518`)

DB
--
* Ensure correct float format when writing to Postgres database (:mr:`3475`)

Structure tools
---------------

* Add atom tagging to ``ase.build.general_surface`` (:mr:`2773`)
* Fix bug where code could return the wrong lattice when trying to fix the handedness of a 2D lattice  (:mr:`3387`)
* Major improvements to :func:`~ase.build.find_optimal_cell_shape`: improve target metric; ensure rotationally invariant results; avoid negative determinants; improved performance via vectorisation (:mr:`3404`, :mr:`3441`, :mr:`3474`). The ``norm`` argument to :func:`~ase.build.supercells.get_deviation_from_optimal_cell_shape` is now deprecated.
* Performance improvements to :class:`ase.spacegroup.spacegroup.Spacegroup` (:mr:`3434`, :mr:`3439`, :mr:`3448`)
* Deprecated :func:`ase.spacegroup.spacegroup.get_spacegroup` as results can be misleading (:mr:`3455`).
  

Calculators / IO
----------------

* Amber: Fix scaling of velocities in restart files (:mr:`3427`)
* Amber: Raise an error if cell is orthorhombic (:mr:`3443`)
* CASTEP

  - **BREAKING** Removed legacy ``read_cell`` and ``write_cell`` functions from ase.io.castep. (:mr:`3435`)
  - .castep file reader bugfix for Windows (:mr:`3379`), testing improved (:mr:`3375`)
  - fix read from Castep geometry optimisation with stress only (:mr:`3445`)

* EAM: Fix calculations with self.form = "eam" (:mr:`3399`)
* FHI-aims
  
  - make free_energy the default energy (:mr:`3406`)
  - add legacy DFPT parser hook (:mr:`3495`)

* FileIOSocketClientLauncher: Fix an unintended API change (:mr:`3453`)
* FiniteDifferenceCalculator: added new calculator which wraps other calculator for finite-difference forces and strains (:mr:`3509`)
* GenericFileIOCalculator fix interaction with SocketIO (:mr:`3381`)
* LAMMPS

  - fixed a bug reading dump file with only one atom (:mr:`3423`)
  - support initial charges (:mr:`2846`, :mr:`3431`)

* MixingCalculator: remove requirement that mixed calculators have common ``implemented_properties`` (:mr:`3480`)
* MOPAC: Improve version-number parsing (:mr:`3483`)
* MorsePotential: Add stress (:mr:`3485`)
* NWChem: fixed reading files from other directories (:mr:`3418`)
* Octopus: Improved IO testing (:mr:`3465`)
* ONETEP calculator: allow ``pseudo_path`` to be set in config (:mr:`3385`)
* Orca: Only parse dipoles if COM is found. (:mr:`3426`)
* Quantum Espresso

  - allow arbitrary k-point lists (:mr:`3339`)
  - support keys from EPW (:mr:`3421`)
  - Fix path handling when running remote calculations from Windows (:mr:`3464`)

* Siesta: support version 5.0 (:mr:`3464`)
* Turbomole: fixed formatting of "density convergence" parameter (:mr:`3412`)
* VASP

  - Fixed a bug handling the ICHAIN tag from VTST (:mr:`3415`)
  - Fixed bugs in CHG file writing (:mr:`3428`) and CHGCAR reading (:mr:`3447`)
  - Fix parsing POSCAR scale-factor line that includes a comment (:mr:`3487`)
  - Support use of unknown INCAR keys (:mr:`3488`)
  - Drop "INCAR created by Atomic Simulation Environment" header (:mr:`3488`)
  - Drop 1-space indentation of INCAR file (:mr:`3488`)
  - Use attached atoms if no atom argument provided to :func:`ase.calculators.vasp.Vasp.calculate` (:mr:`3491`)

GUI
---
* Refactoring of :class:`ase.gui.view.View` to improve API for external projects (:mr:`3419`)
* Force lines to appear black (:mr:`3459`)
* Fix missing Alt+X/Y/Z/1/2/3 shortcuts to set view direction (:mr:`3482`)
* Fix incorrect frame number after using Page-Up/Page-Down controls (:mr:`3481`)
* Fix incorrect double application of `repeat` to `energy` in GUI (:mr:`3492`)

Molecular Dynamics
------------------

* Added Bussi thermostat :class:`ase.md.bussi.Bussi` (:mr:`3350`)
* Added Nose-Hoover chain NVT thermostat :class:`ase.md.nose_hoover_chain.NoseHooverChainNVT` (:mr:`3508`)
* Improve ``force_temperature`` to work with constraints (:mr:`3393`)
* Add ``**kwargs`` to MolecularDynamics, passed to parent Dynamics (:mr:`3403`)
* Support modern Numpy PRNGs in Andersen thermostat (:mr:`3454`)

Optimizers
----------
* **BREAKING** The ``master`` parameter to each Optimizer is now passed via ``**kwargs`` and so becomes keyword-only. (:mr:`3424`)
* Pass ``comm`` to BFGS and CellAwareBFGS as a step towards cleaner parallelism (:mr:`3397`)
* **BREAKING** Removed deprecated ``force_consistent`` option from Optimizer (:mr:`3424`)

Phonons
-------

* Fix scaling of phonon amplitudes (:mr:`3438`)
* Implement atom-projected PDOS, deprecate :func:`ase.phonons.Phonons.dos` in favour of :func:`ase.phonons.Phonons.get_dos` (:mr:`3460`)
* Suppress warnings about imaginary frequencies unless :func:`ase.phonons.Phonons.get_dos` is called with new parameter ``verbose=True`` (:mr:`3461`)

Pourbaix (:mr:`3280`)
---------------------

* New module :mod:`ase.pourbaix` written to replace :class:`ase.phasediagram.Pourbaix`
* Improved energy definition and diagram generation method
* Improved visualisation

Spectrum
--------
* **BREAKING** :class:`ase.spectrum.band_structure.BandStructurePlot`: the ``plot_with_colors()`` has been removed and its features merged into the ``plot()`` method.

Misc
----
* Cleaner bandgap description from :class:`ase.dft.bandgap.GapInfo` (:mr:`3451`)

Documentation
-------------
* The "legacy functionality" section has been removed (:mr:`3386`)
* Other minor improvements and additions (:mr:`2520`, :mr:`3377`, :mr:`3388`, :mr:`3389`, :mr:`3394`, :mr:`3395`, :mr:`3407`, :mr:`3413`, :mr:`3416`, :mr:`3446`, :mr:`3458`, :mr:`3468`)

Testing
-------
* Remove some dangling open files (:mr:`3384`)
* Ensure all test modules are properly packaged (:mr:`3489`)

Units
-----
* Added 2022 CODATA values (:mr:`3450`)
* Fixed value of vacuum magnetic permeability ``_mu0`` in (non-default) CODATA 2018 (:mr:`3486`)

Maintenance and dev-ops
-----------------------
* Set up ruff linter (:mr:`3392`, :mr:`3420`)
* Further linting (:mr:`3396`, :mr:`3425`, :mr:`3430`, :mr:`3433`, :mr:`3469`, :mr:`3520`)
* Refactoring of ``ase.build.bulk`` (:mr:`3390`), ``ase.spacegroup.spacegroup`` (:mr:`3429`)

Earlier releases
================

Releases earlier than ASE 3.24.0 do not have separate release notes and changelog.
Their changes are only listed in the :ref:`releasenotes`.
