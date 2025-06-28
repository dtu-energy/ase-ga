.. _releasenotes:

=============
Release notes
=============

A comprehensive list of changes can be found in the :ref:`changelog`.

Git master branch
=================

:git:`master <>`.

* No changes yet


Version 3.25.0
==============

11 April 2025: :git:`3.25.0 <../3.25.0>`

Version 3.24.0
==============

28 December 2024: :git:`3.24.0 <../3.24.0>`

Requirements
------------

* The minimum supported Python version has increased to 3.9 (:mr:`3473`)

Breaking changes
----------------

* The ``master`` parameter to each Optimizer is now passed via ``**kwargs`` and so becomes keyword-only. (:mr:`3424`)
* Removed legacy ``read_cell`` and ``write_cell`` functions from ase.io.castep. (:mr:`3435`)
* Removed deprecated ``force_consistent`` option from Optimizer (:mr:`3424`)
* :class:`ase.spectrum.band_structure.BandStructurePlot`: the ``plot_with_colors()`` has been removed and its features merged into the ``plot()`` method.

Highlights
----------

* Major improvements to :func:`~ase.build.find_optimal_cell_shape`: improve target metric; ensure rotationally invariant results; avoid negative determinants; improved performance via vectorisation (:mr:`3404`, :mr:`3441`, :mr:`3474`). The ``norm`` argument to :func:`ase.build.supercells.get_deviation_from_optimal_cell_shape` is now deprecated.
* Added new :class:`~ase.calculators.fd.FiniteDifferenceCalculator`, which wraps other calculator for finite-difference forces and strains (:mr:`3509`)
* Added two new MD thermostats: :class:`ase.md.bussi.Bussi` (:mr:`3350`) and :class:`ase.md.nose_hoover_chain.NoseHooverChainNVT` (:mr:`3508`)
* Added atom-projected partial phonon dos to :func:`ase.phonons.Phonons.get_dos` (:mr:`3460`)
* New module :mod:`ase.pourbaix` written to replace
  :class:`ase.phasediagram.Pourbaix` (:mr:`3280`), with improved energy definition and visualisation


Version 3.23.0
==============

31 May 2024: :git:`3.23.0 <../3.23.0>`

* Add :class:`~ase.constraints.FixSubsetCom` to fix the center of mass of the
  specified subset of atoms (:mr:`3193`).

* Add the ``indices`` option to :meth:`ase.Atoms.get_center_of_mass` to compute
  the center of mass of the specified subset of atoms (:mr:`3193`).

* All saddle-point search and minimum-energy path methods have been moved
  into a new :mod:`ase.mep` module. You should start importing from this
  location; e.g., ``from ase.mep import NEB, DyNEB`` (:mr:`2974`).

* Fix :meth:`ase.constraints.FixCom.adjust_momenta` to have no center-of-mass
  momentum (:mr:`2947`)

* :func:`ase.build.surface.graphene` now takes a ``thickness`` argument
  (:mr:`2963`)

* :func:`ase.build.bulk` now assigns initial magnetic moments
  to BCC Fe, Co, and Ni.

* :func:`ase.build.make_supercell` can now control how to order the atoms in
  the supercell via the option ``order`` (:mr:`2800`)

* :meth:`~ase.cell.Cell.mask` returns the mask of nonzero cell vectors,
  an array of three booleans.

* :meth:`~ase.cell.Cell.reciprocal` now guarantees that cell vectors
  which are zero correspond to reciprocal cell vectors that are also
  zero.  Previously the reciprocal cell vectors could have small
  nonzero values due to floating point tolerance.

* The :class:`~ase.cell.Cell` object now has
  :meth:`~ase.cell.Cell.normal` and :meth:`~ase.cell.Cell.normals`
  which calculate normal vectors to one or all pairs of cell vectors.
  Also added
  :meth:`~ase.cell.Cell.area` and
  :meth:`~ase.cell.Cell.areas`, which return the area spanned by one
  or all pairs of cell vectors.

* New ``a2b`` and ``periodic`` formats for :class:`~ase.formula.Formula`
  objects.  See :meth:`ase.formula.Formula.format`.  The ``abc`` format
  has been renamed to ``ab2``.

* IO formats can now be implemented in separate packages and registered
  in ase with the entry point ``ase.ioformats`` in the external package
  configuration. This entry point only accepts objects of the type
  :class:`~ase.utils.plugins.ExternalIOFormat`.

* :class:`ase.phonons.Displacement` now has a ``comm`` keyword argument to
  better support parallelization.

* Fix :func:`ase.build.rotate.minimize_rotation_and_translation` for
  systems with periodic boundary conditions.

* Fix :func:`ase.io.db.row.AtomsRow.charge`, which was looking for
  ``inital_charges`` instead of ``initial_charges``.

* Fixed definition of the Lennard-Jones potential in :mod:`ase.utils.ff`
  (:mr:`2842`)

* The ``x3d`` viewer has improved aesthetics (:mr:`2831`)

* Added a function, :func:`ase.dft.kpoints.mindistance2monkhorstpack`, to
  construct a Monkhorst-Pack grid (:mr:`2811`)

* Fixed turbomole calculator parsing issues when `=` are present in data
  groups (:mr:`2808`)

* Fixed round-trip dict (de)serialization with ``FixedMode`` constraint
  (:mr:`2801`)

* Show total number of rows in ASE database (:mr:`2797`)

* Fix units of magnetic moments in ASE database (:mr:`2793`)

* General improvements to POVRay colors (:mr:`2775`, :mr:`2767`)

* Remove deprecated CLI functions, such as the old ``ase-gui`` command
  (:mr:`2770`)

* Fix calculation of stresses in :func:`ase.stress.full_3x3_to_voigt_6_stress`
  (:mr:`2765`)

* Add polarization parser for FHI-AIMS calculator (:mr:`2763`)

* Fix bug in :func:`ase.atoms.Atoms.get_center_of_mass` (:mr:`2712`)

* Fix multiprocessing bug with :class:`ase.neb.NEB` (:mr:`2695`)

* NEB calculations now call ``get_forces()`` instead of ``get_potential_energy()``
  followed by ``get_forces()`` to prevent duplicate work (:mr:`2678`)

* Improve performance of :func:`ase.build.supercells.make_supercell` for
  large supercells (:mr:`2639`)

* Fix conservation of momentum in :class:`ase.md.langevin.Langevin`
  (:mr:`2632`)

* Fix bug in :class:`ase.vibrations.Vibrations` causing property
  calculations to always use the default method (:mr:`3012`)

* Replaced :class:`ase.phasediagram.Pourbaix` class (to be deprecated)
  with the :mod:`ase.pourbaix` module. The latter includes a
  `~ase.pourbaix.Pourbaix` class able to plot a complete diagram given a
  set of references. The decomposition energy is now shown on a colormap
  and the phase boundaries are determined with a plane intersection method.
  (:mr:`3280`)

Calculators:

* Accelerate :class:`~ase.calculators.emt.EMT` using NumPy more (:mr:`3241`)

* Remove ``ase.calculators.vasp.vasp_auxiliary.xdat2traj``. Use ``ase.io`` instead
  (:mr:`2948`)

* Created new module :mod:`ase.calculators.harmonic` with the
  :class:`ase.calculators.harmonic.HarmonicCalculator`
  for calculations with a Hessian-based harmonic force field. Can be used to
  compute Anharmonic Corrections to the Harmonic Approximation.

* Created new :class:`ase.calculators.plumed.Plumed` that is an interface
  between ASE and Plumed_ for carrying out enhanced sampling methods and MD
  postprocessing.

* :class:`ase.calculators.kim.kimmodel.KIMModelCalculator` updated to allow
  users to change the parameters of OpenKIM portable models at run time (see
  https://openkim.org/doc/repository/kim-content/ for an explanation of types
  of OpenKIM models).

* :class:`ase.calculators.mopac.MOPAC` updated to support MOPAC_
  version 22, which was released under GPL in May 2022.  When using
  this MOPAC version or newer, the output "final heat of formation"
  will be interpreted as potential/free energy for ASE purposes.

* Fix ordering of kpoint and spin indices in
  :func:`ase.calulators.abc.GetOutputsMixin.get_eigenvalues`

* :class:`ase.calculators.vasp.Vasp` now supports the ``efermi`` and ``nelmgw``.

* The :class:`ase.calculators.orca.ORCA` calculator has been
  refactored to be based on :class:`ase.calculators.GenericFileIOCalculator`.

* The :class:`ase.calculators.abinit.AbinitProfile` now correctly executes in
  the user-specified directory.

* Fix reading of eigenvalues in :mod:`ase.calculators.dftb` (:mr:`2819`)

* Added an ``energies`` property to the ``lammpslib`` calculator (:mr:`2799`)

* :func:`ase.calculator.lammpslib.is_upper_triangular` is deprecated. Use the
  method outlined in
  `this StackOverflow post <https://stackoverflow.com/a/26912166>`_ (:mr:`3102`)

* :func:`ase.calculator.lammpslib.convert_cell` is deprecated. Please use
  :func:`ase.calculators.lammps.coordinatetransform.calc_rotated_cell` instead.
  (:mr:`3102`)

* Fix ``lsubsrot`` typo in VASP calculator to ``lsubrot`` (:mr:`2762`)

* Fix atomic polarizability for boron in
  :class:`ase.calculators.bond_polarizability.LippincottStuttman` (:mr:`2759`)

* ``dielectric_tensor`` and ``born_effective_charges`` parsing have been added to
  the aims and VASP calculators (:mr:`2750`)

* Update MOPAC calculator to support v22 (:mr:`2745`)

* Fix ordering of forces in ``lammpsrun`` calculator when used with MPI (:mr:`2737`)

* FLEUR calculator has been removed in favor of external package (:mr:`2736`)

* Compatability fixes for newer versions of Octopus (:mr:`2734`)

* Fixed conversion of force units between LAMMPS and ASE (:mr:`2720`)

* ORCA calculator now returns the correct energy in internal geometry
  optimizations (:mr:`2699`)

* Support for 24Mar2022 version of LAMMPS in ``lammpsrun`` calculator (:mr:`2693`)

* DFTB+ calculator now supports the use of an xTB Hamiltonian (:mr:`2682`)

* Fix WIEN2K writer when ``rmt`` is not ``None`` (:mr:`2677`)

* CP2K calculator support for libcx inputs (:mr:`2676`)

* Fix bug in ``ase run --equation-of-state`` (:mr:`2675`)

* Support PBCs in Plumed calculator (:mr:`2671`)

* Support z-matrix format for ``FixCartesian``` constraints in Siesta (:mr:`2669`)

* Support spin-orbit coupling and non-colinear calculations in Siesta
  (:mr:`2665`, :mr:`2665`)

* Support ``k_grid_density`` keyword in FHI-AIMS calculator (:mr:`2662`)

* Fix storing of ``special_params`` in VASP calculator (:mr:`2661`)

* Fix VASP OUTCAR parsing for systems with over 10 unique elements (:mr:`2658`)

* Cache all available properties in the mixing calculators (:mr:`2646`)

* Refactoring and bug fixes for :class:`ase.calculators.dftd3.DFTD3` (:mr:`2615`)

* Added :func:`ase.calculators.vasp.Vasp.read_vib_modes` and
  :func:`ase.calculators.vasp.Vasp.get_vibrations` to read VASP-calculated
  vibrational modes (:mr:`2605`)

.. _Plumed: https://www.plumed.org/
.. _MOPAC: https://doi.org/10.5281/zenodo.6511958

* :class:`ase.calculators.nwchem.NWChem` now supports
  making initial guesses for wavefunctions using smaller basis
  sets and faster theories

* :class:`ase.calculators.nwchem.NWChem` added support for
  direct and RI-MP2 computations.

Algorithms:

* Fix bugs in :class:`~ase.neighborlist.PrimitiveNeighborList` when used with
  ``bothways=True`` and either ``self_interaction=True`` or ``sorted=True``
  (:mr:`3239`)

* Change :func:`ase.optimize.optimize.Dynamics.irun` so that it yields its
  convergence status *after* logging and dumping for the corresponding step
  (:issue:`1339`; :mr:`3072`)

* :class:`ExpCellFilter` is deprecated. Use :class:`FrechetCellFilter` instead.
  (:mr:`2981`)

* Fix inconsistency of ``mask`` in :class:`FixCartesian` between the argument
  and the attribute. (:mr:`3195`)

Optimizers:

* Add :class:`ase.optimize.climbfixinternals.BFGSClimbFixInternals` class for
  transition state search and optimization along internal reaction coordinates
  (:mr:`2299`)

* :func:`ase.optimize.optimize.Optimizers.irun` and
  :func:`ase.optimize.optimize.Optimizers.run` now respect ``steps=0`` (:issue:`1183`;
  :issue:`1258`; :mr:`2922`).

* Added the ``.trajectory`` attribute to :class:`ase.optimize.optimize.Dynamics`
  (:mr:`2901`).

* Fixed a bug when :class:`ase.optimize.precon.precon.PreconImages` is initialized with
  a list of ``precon`` objects (:mr:`2875`)

* :class:`ase.optimize.mdmin.MDMin` now takes in a ``maxstep`` parameter that
  directly restricts how much atoms can move in a single optimization step
  (:mr:`3065`).

* Removed ``Berny`` optimizer (:mr:`3151`)

Thermochemistry:

* All thermochemistry modules now sort the vibrational energies before
  cutting them down to the physically appropriate amount.

I/O:

* Add support for reading Mulliken, Löwdin, or Hirshfeld atomic charges in
  :func:`ase.io.gaussian.read_gaussian_out` (:mr:`3332`)

* Add support for reading POSCAR files with negative and multiple scaling
  factors in :func:`ase.io.vasp.read_vasp` (:mr:`3029`)

* Fix parsing of dipole moment in :mod:`ase.io.nwchem` from output files
  with trailing spaces (:mr:`3018`)

* Fix to read and write a compressed trajectory file like ``xxx.traj.gz``
  (:mr:`2997`)

* Fix ``index`` of :func:`ase.io.vasp.read_vasp_xdatcar` to be consistent with
  :func:`ase.io.read` (:mr:`2993`)

* Change ``write_lammps_data`` not to reshape the cell by default (:mr:`2986`)

* Fix :func:`ase.io.lammpsdata.write_lammps_data` to tilt the cell correctly
  (:mr:`2986`)

* GaussView-formatted ``.mol`` files can now be read (:mr:`2956`)

* Reading of "chemical json" file types is assigned to names ``*.cjson``
  as used in avogadro2_

* Fixes ``IndexError`` when :func:`lammps_data_to_ase_atoms` is run on a system
  with 1 atom.

* Added several missing ``__init__`` parameters to ``self`` in
  :class:`ase.io.trajectory.TrajectoryReader` and
  :class:`ase.io.trajectory.TrajectoryWriter`.

* Add an option to :func:`ase.io.lammpsdata.write_lamps_data` to print
  the atomic masses.

* Add support for reading CUBE files with "negative" number of atoms, as is common
  in Gaussian.

* Fix parsing of periodic boundary conditions for ``extxyz`` format.

* Increase compatibility of CIF parser

* Extended XYZ writer now works with ``GenericFileIOCalculator``

* Add numerical stress skips for SCF re-initialization in
  :class:`ase.io.aims.AimsOutCalcChunk`.

* Deprecated calling :func:`ase.io.aims.write_aims` with ``velocities``. Use
  ``write_velocities`` instead (:mr:`2910`).

* Fix ``ValueError`` that was raised when using quaternions with
  :func:`ase.io.lammpsrun.lammps_data_to_ase_atoms`.

* :class:`ase.io.castep.CASTEP` now suports reading ``.castep`` force blocks

* :class:`ase.io.castep.CASTEP` no longer reads symmetry operations.

* Added :func:`ase.io.abinit.read_abinit_gsr` to read a netcdf file (:mr:`2855`)

* SDF reader can now read file with >100 atoms (:mr:`2845`)

* :func:`ase.io.v_sim.read_v_sim` now properly accounts for periodic boundary
  conditions (:mr:`2835`)

* Bug fix in the :mod:`ase.io.aims` parser for reading the eigenvalues (:mr:`2825`)

* Recognize ``.castep`` files from Materials Studio (:mr:`2822`)

* Enable reading of GPAW stresses from gpaw-out file (:mr:`2798`)

* Add support for reading a VASP VTST CENTCAR file (:mr:`2751`)

* Fix I/O for FHI-aims when the minimal basis set is used (:mr:`2738`)

* Fix XYZ read/write round-trip when initial and final charges are set (:mr:`2692`)

* CASTEP I/O can now read/write custom masses (:mr:`2686`)

* Default dmol3 extension for arc files is now .arc (:mr:`2686`)

* Add PBC support for :mod:`ase.io.cube` (:mr:`2679`)

* Fix ``sort_by_id=True`` setting in :func:`ase.io.lammpsdata.read_lammps_data`
  (:mr:`2660`)

* More robust GPAW log file reading (:mr:`2609`)

* Fix for PDB file I/O involving trajectories (:mr:`2598`)

* Improved XCrysden file I/O (:mr:`2594`)

* Fix JSON encoder for Atoms objects with ``FixAtoms`` constraints (:mr:`2592`)

* Removed ``ase.io.gaussian_reader`` (:mr:`2329`)

.. _avogadro2: https://www.openchemistry.org/projects/avogadro2


Version 3.22.1
==============

1 December 2021: :git:`3.22.1 <../3.22.1>`

* Fixed compatibility of Brillouin zone plotting with matplotlib 3.5+.


Version 3.22.0
==============

24 June 2021: :git:`3.22.0 <../3.22.0>`


Calculators:

* :class:`ase.calculators.qmmm.ForceQMMM` was updated to enable correct
  handling of various periodic boundary conditions.
  Functions to import and export files with QM/MM mapping were also added.

* It is now possible to use :class:`~ase.calculators.abinit.Abinit`
  together with :class:`~ase.calculators.socketio.SocketIOCalculator`.
  Requires Abinit 9.4+.

* It is now possible to pass a function to
  :class:`~ase.calculators.socketio.SocketIOCalculator` to customize
  startup of a socket client.  This decouples socket I/O calculators
  from :class:`~ase.calculators.calculator.FileIOCalculator`.

* Added :class:`~ase.calculators.socketio.PySocketIOClient`, a helper class
  for using the socket I/O calculator with Python clients.

* OpenKIM calculator updated to support kimpy 2.0.0.

* DFTB+ calculator now reads dipole moments.

Algorithms:

* Dedicated class :class:`ase.vibrations.VibrationsData` to improve
  the representation of vibrational modes and associated data.

* Major refactoring of :class:`ase.vibrations.Vibrations`.
  The calculated vibrational data can now be exported as
  a :class:`~ase.vibrations.VibrationsData` object.

* :meth:`phonons.get_dos` now returns a DOS object based on the new
  framework in :mod:`ase.spectrum`.

* :class:`ase.vibrations.Vibrations` and :class:`ase.phonons.Phonons`
  now use a simplified caching system where forces for each
  displacement are saved in JSON files inside a subdirectory.  This
  breaks old cached calculations.  Old vibrations calculations can be
  ported using a migration tool; see ``python3 -m
  ase.vibrations.pickle2json --help``.

* Added :class:`ase.md.contour_exploration.ContourExploration`.
  It evolves systems at fixed potential energy. This is useful for tracing
  potential energy contour lines or rapidly exploring the potential
  energy surface of a system and can be tuned to preferentially sample
  highly curved regions of the potential energy surface.

* :class:`ase.neb.NEB` has been overhauled and given support for
  preconditioning via a new `precon` argument to its constructor,
  and two newly supported methods, `spline` for spline-interpolated
  tangets and `string` for the string method, both of which support
  preconditioning. The default behaviour should be unchanged.

* Interpolating NEB images on constrained atoms will now raise an
  error if the interpolated positions would become different depending
  on whether the constraints were applied.  Pass
  ``apply_constraint=True`` or ``False`` to
  :meth:`ase.neb.NEB.interpolate` or :func:`ase.neb.interpolate` to
  choose a specific behaviour and silence the error.

* 3D Brillouin zone plots are now guaranteed isometric with Matplotlib 3.3+.

I/O:

* Gaussian input file parsing has been greatly improved.  The parser now
  extracts all variables from the input file.

* Reading of "chemical json" file types with name ``*.cml`` is enabled.

* LAMMPS dump: Reading of elements column added, with priority over types
  if given. All four of the position specifier columns read correctly now.

* Format readers that would by default read or write specific files
  into current working directory no longer do so.  A path, whether
  absolute or relative, is now mandatory for all I/O functions.

* The Siesta .XV format is now a recognized I/O format, ``siesta-xv``.

* Parsing an OUTCAR file will now produce an Atoms object
  with periodic boundary conditions.

Breaking changes:

* For security, ASE no longer uses pickle for any kind of file I/O.
  This is because a maliciously crafted pickle file can execute
  arbitrary code.

  Features that used pickle now either use JSON, no longer support
  saving, or require a manual port of older pickle files using a
  migration tool.  If you have many old calculations and rely on your
  own old (trusted) pickle files which cannot be loaded now, consider
  writing and contributing a migration tool for those files.

  The old PickleTrajectory format can still be loaded
  by manually overriding the security check.

  Pickle is still used for communication between processes started by
  ASE (such as plotting tools in the GUI), which is not a security problem
  since an attacker cannot tamper with the data unless the system is
  already compromised.

GUI:

* Added Finnish translation.

Bug fixes:

* Fix deadlock with DFTD3 calculator in MPI calculations.
* Fix parsing of Quantum Espresso outputs with more than 1000 atoms.
* Write netcdf trajectories compatible with Amber 20.
* Fix bug where constraints could be applied inconsistently in MD
  simulations.
* Allow disabling thermostat and barostat in NPT molecular dynamics.
* Fix problem with whitespace in CIF parser.
* Fix a problem where constraints would be applied inconsistently in
  MD simulations.  As the interactions between MD and constraints are
  not trivial, users should in general verify carefully that simulations
  behave physically correctly.
* Fix issue where occupancies in ``atoms.info`` would subtly change
  type when saved to JSON and reloaded.


Web-page:

* There used to be two versions of the ASE web-page which was quite
  confusing.  The https://wiki.fysik.dtu.dk/ase/dev/ web-page has now been
  dropped.  There is now only https://wiki.fysik.dtu.dk/ase/ and it documents
  the use of the in development version of ASE.


Version 3.21.1
==============

24 January 2021: :git:`3.21.1 <../3.21.1>`

* Fix incorrect positions written to CIF files with mixed boundary
  conditions.
* Writing a CIF with only 1 or 2 lattice vectors will now raise an error since
  CIF cannot represent those systems.
* The name of the Vasp calculator is now ``'vasp'`` as intended.
* Fix attribute error in :meth:`~ase.vibrations.Vibrations.write_jmol`.


Version 3.21.0
==============

18 January 2021: :git:`3.21.0 <../3.21.0>`

General changes:

* :meth:`~ase.Atoms.center` now centers around 0 along directions which
  do not have a cell vector.  Previously this operation had no effect
  in those directions.

* Deprecated the following methods on :class:`~ase.Atoms` as they can
  be replaced by ``~ase.cell.Cell``:
  ``atoms.get_cell_lengths_and_angles()``,
  ``atoms.get_reciprocal_cell()``,
  ``atoms.number_of_lattice_vectors``.
  Instead use ``atoms.cell.cellpar()``, ``atoms.cell.reciprocal()``,
  and ``atoms.cell.rank``, respectively.

* Removed deprecated code on :class:`~ase.Atoms` for handling angles
  in radians.

* :meth:`~ase.Atoms.get_velocities` will now return zeros rather than ``None``
  when there are no velocities, consistently with other optionals such as
  momenta.

* For security reasons, *pickle will no longer be used for persistent
  file storage* in the future.  Pickle has so far been replaced with JSON in
  :class:`ase.io.bundletrajectory.BundleTrajectory`
  and :class:`ase.dft.stm.STM`.
  All remaining use of pickle for persistent storage will be likewise replaced
  in next release.  Users are advised as always not to open pickle-files
  from untrusted sources.

* :func:`ase.utils.opencew` to be replaced by
  :func:`ase.utils.xwopen` which is a contextmanager and ensures
  that the file is closed correctly.

* Clusters created by :mod:`ase.cluster` will no longer have cell vectors
  and will be centered around (0, 0, 0).
  Previously they had a “tight” cell and coordinates
  centered with zero vacuum.

* Refactored external viewers in :mod:`ase.visualize.view`.
  Viewers will now clean up their temporary files correctly on non-UNIX
  platforms.

* Band structure module moved to :mod:`ase.spectrum.band_structure`.

* New objects for working with DOS and collections of DOS in
  :mod:`ase.spectrum`.  To begin with, this will mostly be relevant
  for format readers that want to retrieve such objects from
  calculations.

Command-line interface:

* Added ``ase exec`` sub-command for the :ref:`ase <cli>`
  command line interface.


Algorithms:

* Removed ``ase.build.voids`` (:mr:`2078`)

* Removed unused code in ``ase.transport.tools`` (:mr:`2077`)

* Removed ``ase.visualize.primiplotter`` and ``ase.visualize.fieldplotter``
  (:mr:`2060`)

* Changed units for molecular dynamics modules.  They now accept the
  temperature in Kelvin as a keyword-only argument ``temperature_K``
  and Berendsen NPT accepts the pressure in eV/Å³ as a keyword-only
  argument ``pressure_au``. The previous arguments are still
  available and still take temperature and pressure in whatever unit
  the module used to accept, but now issue a warning.

* Made Andersen thermostat available for molecular dynamics simulation.

* Refactored :class:`ase.neb.NEB`.

* The linear interpolation (:meth:`ase.neb.interpolate`) between images
  now supports cell-interpolation and the use of scaled positions.

* :class:`~ase.neb.SingleCalculatorNEB` is deprecated.  Use
  ``ase.neb.NEB(allow_shared_calculator=True)`` instead.

* Extended constraint :class:`ase.constraints.FixInternals` by
  possibility to fix linear combinations of bond lengths.

* :class:`~ase.constraints.FixInternals` constraints now support
  constraining linear combinations of angles or dihedrals.
  It is also possible to slice atoms objects with
  :class:`~ase.constraints.FixInternals` constraints on them.

* Added :mod:`ase.build.connected` which finds groups of connected
  atoms inside an :class:`~ase.Atoms` object.

* Optimizers and molecular dynamics objects, which may open trajectories
  or logfiles, can now be used as context managers.  Doing so ensures
  correct closing of the files that they open.

* Faster codepath for minimum-image convention (MIC) distance calculations
  with “well-behaved” unit cells.
  This improves the speed of neighbour lists and certain
  constraints.

* Cleanup and deprecations of certain methods on :class:`~ase.phonons.Phonons`.

Calculators:

* Removed ``ase.calculators.ase_qmmm_manyqm`` (:mr:`2092`)

* The ``ignore_bad_restart_file`` argument supported by many calculators
  has been deprecated.  The user should choose this kind of behaviour
  explicitly.

* Cleaned up and fixed multiple issues with
  :class:`~ase.calculators.elk.ELK` calculator.

* Make-shift cleanup and fixes for
  :class:`~ase.calculators.exciting.Exciting` calculator.

* :class:`ase.calculators.abinit.Abinit` updated to work with Abinit 9.

* Improved cleanup of old socket files under some types of failure with
  :class:`ase.calculators.socketio.SocketIOCalculator`.

* :class:`~ase.calculators.vasp.Vasp` now uses the newer implementation
  formerly known as ``Vasp2``. ``Vasp2`` is deprecated.

* Added smooth cutoff option to :class:`ase.calculators.lj.LennardJones`.
  This makes the forces continuous as atoms move past the cutoff radius.

* :class:`~ase.calculators.lj.LennardJones` is now much more efficient.

* Many calculators would change the working directory in order to facilitate
  work with files.  However doing so breaks threading.  This has been fixed
  for most calculators (abinit, lammpsrun, )

I/O:

* Removed ``ase.io.iwm`` (:mr:`2064`)

* Removed ``ase.io.plt`` (:mr:`2057`)

* Reads Wannier90 ``.wout`` files.
  See :func:`ase.io.wannier90.read_wout` and
  :func:`ase.io.wannier90.read_wout_all`.

* :func:`ase.io.pov.write_pov` no longer includes an option to run
  povray on top of the written output.  Instead it returns a renderer
  which can be used like this::

    png_path = write_pov('myfile.pov').render()

* Refactored CIF reader and writer, adding more extensive testing
  and fixing multiple bugs.

* CIF writer now uses up-to-date variable definitions from the CIF
  standard instead of deprecated ones.  Also, it no longer writes columns
  of dummy data that doesn't depend on the atoms.

* Added :class:`ase.io.cif.CIFBlock` for direct access to data
  inside a CIF, and :func:`ase.io.cif.parse_cif` to iterate over
  such blocks from a CIF.

* Fixed many cases of careless I/O handling where format readers or writers
  would open files without necessarily closing them.

* Vasp output formats return atoms with fully periodic boundary conditions
  as appropriate.

* Vasp POSCAR/CONTCAR writer will now use the Vasp5 format by default.

Development:

* Test suite now prints a descriptive header with dependency versions
  including a list of installed/enabled calculators.

* All tests with random numbers now use a specific seed so as to run
  reproducibly.

* CI now supports integration testing with many additional calculators.
  The full list of external calculators that can be integration-tested via CI
  is: Abinit, Asap, CP2K, DFTB, DFTD3, Elk, Espresso,
  Exciting, GPAW, Gromacs, Lammpslib, Lammpsrun, NWChem, Octopus,
  OpenMX, Siesta.


Version 3.20.1
==============

11 August 2020: :git:`3.20.1 <../3.20.1>`

 * Minor fix related to package version requirements on pypi.

Version 3.19.3
==============

11 August 2020: :git:`3.19.3 <../3.19.3>`

 * Minor fix related to package version requirements on pypi.

 * Deprecated calling the :class:`ase.optimize.fire.FIRE` constructor
   with ``maxmove``; please use ``maxstep`` (:mr:`1725`).

Version 3.20.0
==============

8 August 2020: :git:`3.20.0 <../3.20.0>`

General changes:

* Removed old ``ase.data`` modules (:mr:`1720`)

* :meth:`~ase.Atoms.get_calculator` and :meth:`~ase.Atoms.set_calculator`
  are deprecated.  Use ``atoms.calc`` instead.

* ``del atoms.calc`` is deprecated.  Use ``atoms.calc = None`` instead.

* The deprecated ``atoms.cell.pbc`` has been removed.

* More utility methods for ``atoms.symbols``:
  :meth:`~ase.symbols.Symbols.species`,
  :meth:`~ase.symbols.Symbols.search`,
  :meth:`~ase.symbols.Symbols.indices`.

Development:

* Test suite now uses `pytest <https://docs.pytest.org/>`_.
  This means it requires pytest and optionally
  `pytest-xdist <https://github.com/pytest-dev/pytest-xdist>`_ for
  parallelization.  The ``ase test`` command works as before although
  its output will be different and improved.

* Many tests have been improved and simplified, making use of pytest
  for parametrization and test fixtures.

* The continuous integration tests on Gitlab now use custom dockers.
  The docker files can be found at https://gitlab.com/ase/ase-dockers.

* Some calculators can now be tested via Gitlab's CI.

* Code coverage statistics are now available on https://ase.gitlab.io/ase.
  They currently exclude calculators and IO formats.

* Our CI now uses mypy_ for static analysis of the code.

Algorithms:

* Removed ``ase.build.adsorb`` (:mr:`1845`)

* Removed unused code in ``ase.utils.ff`` (:mr:`1844`)

* Removed ``ase.utils.extrapolate`` (:mr:`1808`)

* Functions for attaching structures in :mod:`attach <ase.build>` introduced.

* Standardize optimizers maximum step variable name to maxstep and default
  value to 0.2 for all optimizers.

* Added :class:`Pyberny <ase.optimize.Berny>` geometry optimizer.
  This optimizer can be very
  efficient for molecules, but is currently unreliable and therefore
  should be considered experimental.

* Removed interface to `FindSym
  <https://stokes.byu.edu/iso/findsym.php>`_ due to lack of users and
  maintainers.  If you need this, please find it in git history,
  make it work, and write tests.
  (:mr:`1692`)

* The tangent estimates used to make the nudged elastic band (NEB) plots are
  slightly improved to use center, rather than forward differences. This does
  not affect how NEBs are run; only how they are displayed.

* :class:`ase.neb.NEBTools` now allows the simultaneous plotting of
  all bands from a trajectory of a nudged elastic band calculation (or
  similar); this funciton is also available at the command line as
  ``ase nebplot neb.traj``.

* The image-dependent pair-potential (IDPP) interpolation scheme for
  connecting states---i.e., in a saddle-point search---has been moved
  into the method :func:`ase.neb.idpp_interpolate`. This method is a
  more feature-rich version than that accessible via
  :meth:`ase.neb.NEB.interpolate`.

* Reduced code duplication in the :mod:`ase.ga` module by incorporating the
  'bulk' GA functionality into the corresponding 'standard' modules.
  Using the now deprecated 'bulk' GA modules (i.e.
  :mod:`ase.ga.bulk_startgenerator`, :mod:`ase.ga.bulk_crossovers`,
  :mod:`ase.ga.bulk_mutations` and :mod:`ase.ga.bulk_utilities`) raises
  a warning with pointers to the corresponding 'standard' modules.

* Extended the genetic algorithm to cases where 1 or 2 cell vectors are
  part of the global optimization problem, which can be useful in searching
  for nanowire and thin film structures.

* Added a new tutorial on molecular crystal structure prediction using
  a genetic algorithm, see :ref:`ga_molecular_crystal_tutorial`.

* Allow setting the initial hessian in `optimize.BFGS` via the keyword `alpha` or
  explicitly via `opt.H0 = ...` after instantiation.

Command-line interface:

* New ``dimensionality`` sub-command for the :ref:`ase <cli>` command line
  interface.

* Added a diff CLI for displaying and comparing the positions, forces,
  and energies of atoms objects.  The classes and functions used in
  the CLI are also available in :mod:`ase.cli.template`.

I/O:

* The ``ase db db1.db <selection> --insert-into db2.db`` command now respects
  ``--limit`` and ``--offset``.

* Read and write support for qball sys file format.

* Added write support for the Vasp 5 XDATCAR file format.

* Added Z-matrix parser for use in input/output file readers.

* Added support for writing prismatic and computem xyz file. Required arguments
  to write mustem xtl file have been updated to be consistent with prismatic
  and computem xyz file export.

* Removed ETSF format reader since it depends on ScientificPython
  which requires Python 2.7.

* Removed Dacapo-NetCDF reader which has not worked since ancient times.
  (:mr:`1892`)

GUI:

* Use Ctrl+C, Ctrl+X, and Ctrl+V to copy/cut/paste atoms using the
  operating system's clipboard.  The copies use the JSON format.

* Removed old GUI modules which were never fully ported to Tkinter.
  If you miss them, please find them in git history and rehabilitate
  them.

 * Added Russian translation.

Calculators:

* Multiple improvements and bugfixes to OpenMX calculator;
  OpenMX calculator now supports OpenMX 3.9.

* Added :mod:`ORCA <ase.calculators.orca>` calculator.

* Added :mod:`GAMESS-US <ase.calculators.gamess_us>` calculator.

* Removed interface to :ref:`Dacapo <jacapo>` due to lack of users and
  maintainers.
  (:mr:`1721`, :mr:`1604`)

* Completely refactored :mod:`Gaussian <ase.calculators.gaussian>` calculator.
  The new calculator should be completely backwards compatible with the
  previous one, while having a more flexible design and supporting more
  keyword arguments.

* Added :mod:`GaussianOptimizer <ase.calculators.gaussian>` and
  :mod:`GaussianIRC <ase.calculators.gaussian>` classes for performing geometry
  optimization and IRC calculations with the Gaussian calculator. These
  classes are the canonical way to use Gaussian's built-in geometry
  optimization routines.

* Fixed ``kpts`` option of :class:`ase.calculators.espresso.Espresso`
  so that specifying a Γ-point calculation with ``kpts=(1, 1, 1)``
  does not enable the optimized codepath (which halves memory and
  cpu). Use ``kpts=None`` to enable the optimized codepath.

* Added the properties ``stresses`` and ``energies`` to the
  Lennard-Jones potential :class:`ase.calculators.lj.LennardJones`.
  Functionality for other properties should be unchanged. Testing and
  documentation have also been expanded.


.. _mypy: http://mypy-lang.org/


Version 3.19.2
==============

22 July 2020: :git:`3.19.2 <../3.19.2>`

* Compatibility fixes related to matplotlib:
  Update png writer to be compatible with matplotlib 3.3.0.
  Update incompatible calls to ``matplotlib.use()``.

* Deprecated calling :class:`ase.calculators.vasp.vasp2.Vasp2` constructor
  with directory in ``label`` parameter (:mr:`1940`)

Version 3.19.1
==============

4 April 2020: :git:`3.19.1 <../3.19.1>`

* Update png writer to be compatible with matplotlib 3.2.


Version 3.19.0
==============

16 December 2019: :git:`3.19.0 <../3.19.0>`

General changes:

* :func:`ase.build.bulk` now supports elements with tetragonal and
  rhombohedral lattices.

* The ``rank`` and ``size`` constants from the :mod:`ase.parallel` module have
  been deprecated.  Use ``world.rank`` and ``world.size`` instead
  (and ``from ase.parallel import world``).

* ``atoms.set_masses('most_common')`` now sets the masses of each
  element according to most common isotope as stored in
  ``ase.data.atomic_masses_common``.

* :mod:`ase.utils.parsemath` added to utils. This module parses simple
  mathematical expressions and returns their numerical value.

* Plotting functions (such as band structure, EOS, ...)
  no longer show the figure by default.

* :class:`~ase.Atoms` constructor now accepts ``velocities`` as keyword.

* Documentation: New set of :ref:`introductory ASE tutorials <gettingstarted>`.

* More detailed output of ``ase info --formats``.

* For completeness, :mod:`ase.lattice` now also supports the 1D
  Bravais lattice.

Algorithms:

* Added :class:`~ase.md.analysis.DiffusionCoefficient` so one can
  calculate atom/molecule mobility from trajectory as a function of
  time.

* Added general linear parametric constraints :class:`ase.constraints.FixParametricRelations`,
  :class:`ase.constraints.FixScaledParametricRelations`, and
  :class:`ase.constraints.FixCartesianParametricRelations` to
  :mod:`ase.constraints`. These constraints are based off the work
  in: :arxiv:`1908.01610`, and allows for the positions and cell of a
  structure to be optimized in a reduced parameter space.

* Added :func:`ase.build.graphene` for building graphene monolayers.

* Added :mod:`ase.md.switch_langevin` module for thermodynamic
  integration via MD simulations.

* Implemented "dynamic" or "ideal gas" contribution from atomic
  momenta to stress tensor Use :meth:`<ase.Atoms.get_stress>`, e.g.,
  ``atoms.get_stress(include_ideal_gas=True)``.

Calculators:

* Added :mod:`Q-Chem <ase.calculators.qchem>` calculator.

* Added :class:`~ase.calculators.psi4.Psi4` calculator.

* Added :class:`~ase.calculators.demonnano.DemonNano` calculator.

* Added :mod:`OpenKIM <ase.calculators.kim>` calculator,
  a special calculator for `OpenKim <https://openkim.org/>`_ models.

* Gulp calculator now provides stress tensor.

* The :mod:`NWChem <ase.calculators.nwchem>` calculator has been completely rewritten, and now supports
  `DFT <https://github.com/nwchemgit/nwchem/wiki/Density-Functional-Theory-for-Molecules>`_,
  `SCF (Hartree Fock) <https://github.com/nwchemgit/nwchem/wiki/Hartree-Fock-Theory-for-Molecules>`_,
  `MP2 <https://github.com/nwchemgit/nwchem/wiki/MP2>`_,
  `CCSD <https://github.com/nwchemgit/nwchem/wiki/CCSD>`_,
  and `TCE <https://github.com/nwchemgit/nwchem/wiki/TCE>`_ calculations with gaussian-type orbitals.
  The calculator also now supports
  `plane-wave calculations <https://github.com/nwchemgit/nwchem/wiki/Plane-Wave-Density-Functional-Theory>`_,
  including band structure calculations through ASE's :class:`~ase.dft.band_structure.BandStructure` utilities.
  To facilitate these changes, the format of the calculator keywords has been changed. Please read the updated
  :mod:`NWChem <ase.calculators.nwchem>` calculator documentation for more details.

* :class:`~ase.calculators.siesta.siesta.Siesta` calculator refactored.
  The Siesta calculator now supports the band structure machinery.
  There is only a single Siesta calculator now covering all versions of Siesta,
  consistently with other ASE calculators.

* Added :mod:`~ase.calculators.mixing` module for the linear
  combination of arbitrary :mod:`~ase.calculators`.

* New :class:`ase.calculators.idealgas.IdealGas` calculator for
  non-interacting atoms.  The calculator does nothing.  This can be
  useful for testing.

* :class:`~ase.calculators.emt.EMT` calculator now support
  atom-specific energies as per ``atoms.get_energies()``.

I/O:

* Read and write support for RMCProfile (rmc6f) file format.

* Write support for Materials Studio xtd files.

* More efficient storage of the "data" part of rows in the :mod:`ase.db`
  database.  NumPy arrays are now stored in binary format instead of as text
  thereby using approximately a factor of two less space when storing numbers
  of ``np.float64``.

* The :mod:`~ase.io.pov` module can now render high-order bonds.

* :class:`~ase.Atoms` now provides the general-purpose JSON mechanism
  from :mod:`ase.io.jsonio`.

* Added :mod:`ase.data.pubchem` module to search for structures
  in the `PubChem <https://pubchem.ncbi.nlm.nih.gov/>`_ database.

GUI:

* It is now possible to copy and paste atoms: The "add atoms" function
  (Ctrl+A) will suggest the atoms in the current selection by default.


Version 3.18.2
==============

15 December 2019: :git:`3.18.2 <../3.18.2>`

* Fix an issue with the binary package (wheel) of 3.18.1.
  No bugfixes as such.

* Deprecated ``ase.calculator.siesta.base_siesta``


Version 3.18.1
==============

20 September 2019: :git:`3.18.1 <../3.18.1>`

* Multiple bugfixes.  Most importantly, deprecate ``atoms.cell.pbc``
  in order to avoid complexities from dealing with two
  ways of manipulating this piece of information.
  Use ``atoms.pbc`` instead; this works the same as always.
  Also, the :class:`~ase.cell.Cell` object now exposes almost the entire
  ``ndarray`` interface.  For a list of smaller bugfixes, see the git log.

* Deprecated ``ase.Atoms.get_number_of_atoms`` (:mr:`1295`)


Version 3.18.0
==============

19 July 2019: :git:`3.18.0 <../3.18.0>`

General changes:

* ASE no longer supports Python2.

* ``atoms.cell`` is now a :class:`~ase.cell.Cell` object.
  This object resembles a 3x3 array and also provides shortcuts to many common
  operations.

* Preliminary :class:`~ase.formula.Formula` type added.  Collects all
  formula manipulation functionality in one place.

* :class:`~ase.symbols.Symbols` objects, like ``atoms.symbols``, now have a
  :attr:`~ase.symbols.Symbols.formula` attribute.

* Added classes to represent primitive Bravais lattices and data
  relating to Brillouin zones to :mod:`ase.lattice`.  Includes 2D
  lattices.

* New :class:`~ase.dft.kpoints.BandPath` class to represent a band path
  specification like ``'GXL'`` along with actual k-point coordinates.
  :class:`~ase.dft.band_structure.BandStructure` objects now have a band
  path.

* :func:`ase.dft.kpoints.bandpath` now returns a
  :class:`~ase.dft.kpoints.BandPath` object.  Generation
  of band paths now works for (almost) any cell.

* Use ``atoms.cell.bandpath()`` as a shortcut to generate band paths.

* New holonomic :class:`constraint <ase.constraints.FixLinearTriatomic>`
  for trilinear molecules.

* Added ``ase info --calculators`` option which shows a list of
  calculators and whether they appear to be installed.

* Added :func:`ase.build.surfaces_with_termination.surfaces_with_termination`,
  a tool to build surfaces with a particular termination.

* Use the shortcut ``with ase.utils.workdir('mydir', mkdir=True):
  <code>`` to temporarily change directories.

* The ``ase test`` command now properly autocompletes test names and
  calculator names.

* Added keyword, ``atoms.wrap(pretty_translation=True)``, to minimize
  the scaled positions of the atoms.

Calculators:

* Added interface to :mod:`ACE-Molecule <ase.calculators.acemolecule>`.

* NWChem calculator now supports TDDFT runs.

* Multiple improvements to the ONETEP Calculator. Input files can now be
  written that specify LDOS, bsunfolding and many other functionalities.

* Calculation of stress tensor implemented for
  :class:`~ase.calculators.emt.EMT` potential.

* The :class:`~ase.calculators.octopus.Octopus` calculator now
  provides the stress tensor.

* Reworked :class:`~ase.calculators.lammpsrun.LAMMPS` calculator.  The
  calculator should now behave more consistently with other ASE
  calculators.

* Gromacs calculator updated to work with newer Gromacs.

* Fleur calculator updated to work with newer Fleur.

* Added :class:`~ase.calculators.ACN`, a QM/MM forcefield for acetonitrile.

* Improved eigenvalue parsing with Siesta calculator.

Algorithms:

* Determine Bravais lattice for any 2D or 3D cell using
  ``atoms.cell.get_bravais_lattice()``.

* Added function to Minkowski reduce a cell.

* Improved stability of Niggli reduction algorithm.

* Supercell generation using ``ase.build.make_supercell()`` now uses
  a constructive algorithm instead of cutting which was prone to tolerance
  errors.

* Setting an MD velocity distribution now preserves the temperature
  by default.

* :class:`Analysis tool <ase.geometry.analysis.Analysis>` for extracting
  bond lengths and angles from atoms.

* Dynamics and structure optimizers can now run as an iterator using the
  new ``irun()`` mechanism::

    for conv in opt.irun(fmax=0.05):
        print('hello')

  This makes it easier to execute custom code during runs.  The ``conv``
  variable indicates whether the current iteration meets the convergence
  criterion, although this behaviour may change in future versions.

* The genetic algorithm module :mod:`ase.ga` now has operators for crystal
  structure prediction. See :ref:`ga_bulk_tutorial`.

* New :func:`ase.geometry.dimensionality.analyze_dimensionality` function.
  See: :ref:`dimtutorial`.

* New :func:`ase.utils.deltacodesdft.delta` function:  Calculates the
  difference between two DFT equation-of-states.  See the new :ref:`dcdft tut`
  tutorial.

* Holonomic :class:`~ase.constraints.FixLinearTriatomic` for QM/MM
  calculations.

* The :class:`~ase.neighborlist.NeighborList` now uses kdtree from Scipy
  for improved performance.  It also uses Minkowsky reduction
  to improve performance for unusually shaped cells.

I/O:

* Database supports user defined tables

* Support for reading and writing DL_POLY format.

* Support for reading CP2K DCD format.

* Support for EON .con files with multiple images.

* Support for writing Materials Studio xtd format.

* Improved JSON support.  :ref:`cli` tools like :program:`ase
  band-structure` and :program:`ase reciprocal` now work with
  JSON representations of band structures and paths.

* Support reading CIF files through the
  `Pycodcif <http://wiki.crystallography.net/cod-tools/CIF-parser/>`_
  library.  This can be useful for CIF features that are not supported
  by the internal CIF parser.

* :ref:`MySQL and MariaDB <MySQL_server>` are supported as database backend

* Support for writing isosurface information to POV format
  with :func:`ase.io.pov.add_isosurface_to_pov`

GUI:

 * Quickinfo dialog automatically updates when switching image.

 * Display information about custom arrays on Atoms objects; allow colouring
   by custom arrays.

 * Improved color scales.


Version 3.17.0
==============

12 November 2018: :git:`3.17.0 <../3.17.0>`

General changes:

* ``atoms.symbols`` is now an array-like object which works
  like a view of ``atoms.numbers``, but based on chemical symbols.
  This enables convenient shortcuts such as
  ``mask = atoms.symbols == 'Au'`` or
  ``atoms.symbols[4:8] = 'Mo'``.

* Test suite now runs in parallel.

* New :class:`~ase.dft.pdos.DOS` object for representing and plotting
  densities of states.

* Neighbor lists can now :meth:`get connectivity matrices
  <ase.neighborlist.NeighborList.get_connectivity_matrix>`.

* :ref:`ase convert <cli>` now provides options to execute custom code
  on each processed image.

* :class:`~ase.phonons.Phonons` class now uses
  the :class:`~ase.dft.pdos.DOS` and
  :class:`~ase.dft.band_structure.BandStructure` machinery.

* Positions and velocities can now be initialized from phononic
  force constant matrix; see
  :func:`~ase.md.velocitydistribution.PhononHarmonics`.

Algorithms:

* New Gaussian Process (GP) regression optimizer
  (:class:`~ase.optimize.GPMin`).  Check out this `performance test
  <https://gpaw.readthedocs.io/devel/ase_optimize/ase_optimize.html>`_.

* New filter for lattice optimization,
  :class:`~ase.constraints.ExpCellFilter`, based on an exponential
  reformulation of the degrees of freedom pertaining to the cell.
  This is probably significantly faster than
  :class:`~ase.constraints.UnitCellFilter`.

* :class:`~ase.constraints.UnitCellFilter` now supports scalar pressure and
  hydrostatic strain.

* Compare if two bulk structure are symmetrically equivalent with
  :class:`~ase.utils.structure_comparator.SymmetryEquivalenceCheck`.

* :class:`~ase.neb.NEB` now supports a boolean keyword,
  ``dynamic_relaxation``, which will freeze or unfreeze images
  according to the size of the spring forces so as to save
  force evaluations.  Only implemented for serial NEB calculations.

* Writing a trajectory file from a parallelized :class:`~ase.neb.NEB`
  calculation is now much simpler.  Works the same way as for the serial
  case.

* New :class:`~ase.constraints.FixCom` constraint for fixing
  center of mass.

Calculators:

* Added :class:`ase.calculators.qmmm.ForceQMMM` force-based QM/MM calculator.

* Socked-based interface to certain calculators through the
  :mod:`~ase.calculators.socketio` module:
  Added support for
  communicating coordinates, forces and other quantities over
  sockets using the i-PI protocol.  This removes the overhead for
  starting and stopping calculators for each geometry step.
  The calculators which best support this feature are Espresso,
  Siesta, and Aims.

* Added calculator for :mod:`OpenMX <ase.calculators.openmx>`.

* Updated the :class:`~ase.calculators.castep.Castep` calculator as well as
  the related I/O methods in order to be more forgiving and less reliant on
  the presence of a CASTEP binary. The ``castep_keywords.py`` file has been
  replaced by a JSON file, and if its generation fails CASTEP files can still
  be read and written if higher tolerance levels are set for the functions that
  manipulate them.

* :class:`~ase.calculators.espresso.Espresso`
  and :mod:`~ase.calculators.dftb` now support the
  :class:`~ase.dft.band_structure.BandStructure` machinery
  including improved handling of kpoints, ``get_eigenvalues()``,
  and friends.

I/O:

* CIF reader now parses fractional occupancies if present.
  The GUI visualizes fractional occupancies in the style of Pacman.

* Support for downloading calculations from the Nomad archive.
  Use ``ase nomad-get nmd://<uri> ...`` to download one or more URIs
  as JSON files.  Use the :mod:`ase.nomad` module to download
  and work with Nomad entries programmatically.  ``nomad-json``
  is now a recognized IO format.

* Sequences of atoms objects can now be saved as animations using
  the mechanisms offered by matplotlib.  ``gif`` and ``mp4`` are now
  recognized output formats.

Database:

* The :meth:`ase.db.core.Database.write` method now takes a ``id`` that
  allows you to overwrite an existing row.

* The :meth:`ase.db.core.Database.update` can now update the Atoms and the data
  parts of a row.

* The :meth:`ase.db.core.Database.update` method will no longer accept a list of
  row ID's as the first argument.  Replace this::

      db.update(ids, ...)

  with::

      with db:
          for id in ids:
              db.update(id, ...)

* New ``--show-keys`` and ``--show-values=...`` options for the
  :ref:`ase db <cli>` command line interface.

* Optimized performance of ase db, with enhanced speed of
  queries on key value pairs for large SQLite (.db) database files.
  Also, The ase db server (PostgreSQL) backend now uses
  native ARRAY and JSONB data types for storing NumPy arrays and
  dictionaries instead of the BYTEA datatype. Note that backwards
  compatibility is lost for the postgreSQL backend, and that
  postgres version 9.4+ is required.

GUI:

* Added callback method :meth:`ase.gui.gui.GUI.repeat_poll` to the GUI.
  Useful for programmatically updating the GUI.

* Improved error handling and communication with subprocesses (for plots)
  in GUI.

* Added Basque translation.

* Added French translation.

Version 3.16.2
==============

4 June 2018: :git:`3.16.2 <../3.16.2>`

* Fix test failure for newer versions of flask due to error within the test itself.  Fix trajectory format on bigendian architectures.  Fix issue with trajectory files opened in append mode where header would not be written correctly for images with different length, atomic species, boundary conditions, or constraints.


Version 3.16.0
==============

21 March 2018: :git:`3.16.0 <../3.16.0>`

* New linear-scaling neighbor list
  available as a function :meth:`~ase.neighborlist.neighbor_list`.

* Castep calculator: option for automatic detection of pseudopotential files from a given directory (castep_pp_path); support for GBRV pseudopotential library; updated outfile parsing to comply with CASTEP 18.1.

* New LAMMPS calculator LAMMPSlib utilizing the Python bindings provided by LAMMPS instead of file I/O. Very basic calculator but can serve as base class for more sophisticated ones.

* Support for µSTEM xtl data format.

* New scanning tunnelling spectroscopy (STS) mode for
  :class:`~ase.dft.stm.STM` simulations.

* New method, :meth:`~ase.Atoms.get_angles`, for calculating multiple angles.

* New ``ase reciprocal`` :ref:`command <cli>` for showing the
  1. Brilluin zone, **k**-points and special points.

* New ``ase convert`` :ref:`command <cli>` for converting between file formats.

* Improved XRD/SAXS module:  :mod:`ase.utils.xrdebye`.

* New cell editor for the GUI.

* Improved "quick info" dialog in the GUI.  The dialog now lists results
  cached by the calculator.

* The "add atoms" dialog now offers a load file dialog as was the case before the tkinter port.  It also provides a chooser for the G2 dataset.

* Interface for the :mod:`CRYSTAL <ase.calculators.crystal` code has been
  added.

* The :func:`ase.dft.bandgap.bandgap` function used with ``direct=True``
  will now also consider spin-flip transitions.  To get the spin-preserving
  direct gap (the old behavior), use::

      min(bandgap(..., spin=s, direct=True) for s in [0, 1])

* Bug fixed in the :meth:`ase.phonons.Phonons.symmetrize` method when using an
  even number of repeats.


Version 3.15.0
==============

28 September 2017: :git:`3.15.0 <../3.15.0>`

* If you are running your Python script in :mod:`parallel <ase.parallel>`
  then by default, :func:`ase.io.read` and :func:`ase.io.iread` will read on
  the master and broadcast to slaves, and :func:`ase.io.write` will only
  write from master.  Use the new keyword ``parallel=False`` to read/write
  from the individual slaves.

* New ``ase find`` :ref:`command <cli>` for finding atoms in files.

* Added :class:`Espresso <ase.calculators.espresso.Espresso>` calculator for
  Quantum ESPRESSO in module :mod:`ase.calculators.espresso`.

* The :func:`ase.dft.kpoints.get_special_points` function has a new call
  signature:  Before it was ``get_special_points(lattice, cell)``, now it is
  ``get_special_points(cell, lattice=None)``.  The old way still works, but
  you will get a warning.

* The :class:`ase.dft.dos.DOS` object will now use linear tetrahedron
  interpolation of the band-structure if you set ``width=0.0``.  It's slow,
  but sometimes worth waiting for.  It uses the
  :func:`ase.dft.dos.linear_tetrahedron_integration` helper function.

* :func:`ase.io.read` can now read QBox output files.

* The :mod:`ase.calculators.qmmm` module can now also use
  :ref:`Turbomole <turbomole qmmm>` and :mod:`DFTB+ <ase.calculators.dftb>`
  as the QM part.

* New :ref:`db tutorial` tutorial.

* :mod:`ase.gui`:  Improved atom colouring options; support the Render Scene (povray) and Ctrl+R rotation features again; updated German and Chinese translations.

* Get the :class:`~ase.spacegroup.Spacegroup` object from an
  :class:`~ase.Atoms` object with the new :func:`ase.spacegroup.get_spacegroup`
  function.


Version 3.14.1
==============

28 June 2017: :git:`3.14.1 <../3.14.1>`.

* Calling the :func:`ase.dft.bandgap.bandgap` function with ``direct=True``
  would return band indices that were off by one.  Fixed now.


Version 3.14.0
==============

20 June 2017: :git:`3.14.0 <../3.14.0>`.

* Python 2.6 no longer supported.

* The command-line tools :program:`ase-???` have been replaced by a
  single :program:`ase` command with sub-commands (see :ref:`cli`).
  For help, type::

      $ ase --help
      $ ase sub-command --help

* The old :program:`ase-build` command which is now called
  :program:`ase build` will no longer add vacuum by default.  Use
  ``ase build -V 3.0`` to get the old behavior.

* All methods of the :class:`~ase.Atoms` object that deal with angles now
  have new API's that use degrees instead of radians as the unit of angle
  (:meth:`~ase.Atoms.get_angle`, :meth:`~ase.Atoms.set_angle`,
  :meth:`~ase.Atoms.get_dihedral`, :meth:`~ase.Atoms.set_dihedral`,
  :meth:`~ase.Atoms.rotate_dihedral`, :meth:`~ase.Atoms.rotate`,
  :meth:`~ase.Atoms.euler_rotate`).

  The old way of calling these methods works as always, but will give
  you a warning.  Example:

  >>> water.get_angle(0, 1, 2)  # new API
  104.52
  >>> water.get_angle([0, 1, 2])  # old API
  /home/jensj/ase/ase/atoms.py:1484: UserWarning: Please use new API (which will return the angle in degrees): atoms_obj.get_angle(a1,a2,a3)*pi/180 instead of atoms_obj.get_angle([a1,a2,a3])
  1.8242181341844732

  Here are the changes you need to make in order to get rid of warnings:

  Old API:

  >>> a1 = atoms.get_angle([0, 1, 2])
  >>> atoms.set_angle([0, 1, 2], pi / 2)
  >>> a2 = atoms.get_dihedral([0, 1, 2, 3])
  >>> atoms.set_dihedral([0, 1, 2, 3], pi / 6)
  >>> atoms.rotate_dihedral([0, 1, 2, 3], 10.5 * pi / 180)
  >>> atoms.rotate('z', pi / 4)
  >>> atoms.rotate_euler(phi=phi, theta=theta, psi=psi)

  New API:

  >>> a1 = atoms.get_angle(0, 1, 2) * pi / 180
  >>> atoms.set_angle(0, 1, 2, angle=90)
  >>> a2 = atoms.get_dihedral(0, 1, 2, 3) * pi / 180
  >>> atoms.set_dihedral(0, 1, 2, 3, angle=30)
  >>> atoms.rotate_dihedral(0, 1, 2, 3, angle=10.5)
  >>> atoms.rotate(45, 'z')
  >>> atoms.euler_rotate(phi=phi * 180 / pi,
  ...                    theta=theta * 180 / pi,
  ...                    psi=psi * 180 / pi)

* The web-interface to the :mod:`ase.db` module now uses Bootstrap and looks
  much nicer.  Querying the database is also much easier.  See
  https://cmrdb.fysik.dtu.dk for an example.

* The PostgreSQL backend for :mod:`ase.db` can now contain more than one ASE
  database.

* An ASE database can now have :ref:`metadata` describing the data.
  Metadata is a dict with any of the following keys: ``title``,
  ``key_descriptions``, ``default_columns``, ``special_keys`` and
  ``layout``.

* :data:`ase.data.atomic_masses` has been updated to IUPAC values from
  2016. Several elements will now have different weights which will affect
  dynamic calculations. The old values can be recovered like this:

  >>> from ase.data import atomic_masses_legacy
  >>> atoms.set_masses(atomic_masses_legacy[atoms.numbers])

* New :func:`ase.data.isotopes.download_isotope_data` function for getting
  individual isotope masses from NIST.

* New :func:`ase.eos.calculate_eos` helper function added.

* Added DeltaCodesDFT data: :data:`ase.collections.dcdft`.

* :mod:`ase.gui` can now load and display any sequence of :class:`~ase.Atoms`
  objects; it is no longer restricted to sequences with a constant number
  of atoms or same chemical composition.

* Trajectory files can now store any sequence of :class:`~ase.Atoms`
  objects.  Previously, atomic numbers, masses, and constraints were
  only saved for the first image, and had to apply for all subsequent ones.

* Added calculator interface for DMol\ :sup:`3`.

* Added calculator interface for GULP.

* Added file formats .car, .incoor, and .arc, related to DMol\ :sup:`3`.

* New function for interpolating from Monkhors-Pack sampled values in the BZ
  to arbitrary points in the BZ:
  :func:`ase.dft.kpoints.monkhorst_pack_interpolate`.

* New *band-structure* command for the :program:`ase` :ref:`cli`.

* Two new functions for producing chemical formulas:
  :func:`ase.utils.formula_hill` and :func:`ase.utils.formula_metal`.

* The :func:`ase.dft.bandgap.get_band_gap` function is now deprecated.  Use
  the new one called :func:`ase.dft.bandgap.bandgap` (it's more flexible and
  returns also band indices).

* New :mod:`Viewer for Jupyter notebooks <ase.visualize.nglview>`.


Version 3.13.0
==============

7 February 2017: :git:`3.13.0 <../3.13.0>`.

* The default unit-cell when you create an :class:`~ase.Atoms` object has
  been changed from ``[[1,0,0],[0,1,0],[0,0,1]]`` to
  ``[[0,0,0],[0,0,0],[0,0,0]]``.

* New :attr:`ase.Atoms.number_of_lattice_vectors` attribute equal to,
  big surprise, the number of non-zero lattice vectors.

* The :meth:`ase.Atoms.get_cell` method has a new keyword argument
  ``complete``.  Use ``atoms.get_cell(complete=True)`` to get a complete
  unit cell with missing lattice vectors added at right angles to the
  existing ones.  There is also a function :func:`ase.geometry.complete_cell`
  that will complete a unit cell.

* :func:`~ase.build.graphene_nanoribbon` no longer adds 2.5 Å of vacuum by
  default.

* All functions that create molecules, chains or surfaces
  (see the :mod:`ase.build` module) will no longer add "dummy" lattice
  vectors along the non-periodic directions.  As an example, the surface
  functions will generate unit cells of the type
  ``[[a1,a2,0],[b1,b2,0],[0,0,0]]``.  In order to define all three lattice
  vectors, use the ``vacuum`` keyword that all
  of the 0-d, 1-d and 2-d functions have or, equivalently, call the
  :meth:`~ase.Atoms.center` method.

* Many of the :ref:`surface generating functions <surfaces>` have changed
  their behavior when called with ``vacuum=None`` (the default).  Before, a
  vacuum layer equal to the interlayer spacing would be added on the upper
  surface of the slab. Now, the third axis perpendicular to the surface will be
  undefined (``[0, 0, 0]``).  Use ``vacuum=<half-the-interlater-distance>`` to
  get something similar to the old behavior.

* New :func:`ase.geometry.is_orthorhombic` and
  :func:`ase.geometry.orthorhombic` functions added.

* :mod:`ase.gui` now works on Python 3.

* NEB-tools class has been renamed to :class:`~ase.neb.NEBTools`.

* :mod:`Optimizers <ase.optimize>` now try force-consistent energies if
  possible (instead of energies extrapolated to 0.0 K).


Version 3.12.0
==============

24 October 2016: :git:`3.12.0 <../3.12.0>`.

* New :class:`ase.constraints.ExternalForce` constraint.

* Updated :mod:`ase.units` definition to CODATA 2014. Additionally, support
  for older versions of CODATA was added such that the respective units can
  be created by the user when needed (e.g. interfacing codes with different
  CODATA versions in use).

* New :mod:`ase.calculators.checkpoint` module.  Adds restart and rollback
  capabilities to ASE scripts.

* Two new flawors of :class:`~ase.neb.NEB` calculations have been added:
  ``method='eb'`` and ``method='improvedtangent'``.

* :func:`ase.io.write` can now write XSD files.

* Interfaces for deMon, Amber and ONETEP added.

* New :ref:`defects` tutorial and new super-cell functions:
  :func:`~ase.build.get_deviation_from_optimal_cell_shape`,
  :func:`~ase.build.find_optimal_cell_shape`,
  :func:`~ase.build.make_supercell`.

* New :class:`~ase.dft.band_structure.BandStructure` object.  Can identify
  special points and create nice plots.

* Calculators that inherit from :class:`ase.calculators.calculator.Calculator`
  will now have a :meth:`~ase.calculators.calculator.Calculator.band_structure`
  method that creates a :class:`~ase.dft.band_structure.BandStructure` object.

* Addition to :mod:`~ase.geometry` module:
  :func:`~ase.geometry.crystal_structure_from_cell`.

* New functions in :mod:`ase.dft.kpoints` module:
  :func:`~ase.dft.kpoints.parse_path_string`,
  :func:`~ase.dft.kpoints.labels_from_kpts` and
  :func:`~ase.dft.kpoints.bandpath`.

* Helper function for generation of Monkhorst-Pack samplings and BZ-paths:
  :func:`ase.calculators.calculator.kpts2ndarray`.

* Useful class for testing band-structure stuff:
  :class:`ase.calculators.test.FreeElectrons`.

* The ``cell`` attribute of an :class:`~ase.Atoms` object and the ``cell``
  keyword for the :class:`~ase.Atoms` constructor and the
  :meth:`~ase.Atoms.set_cell` method now accepts unit cells given ase
  ``[a, b, c, alpha, beta, gamma]``, where the three angles are in degrees.
  There is also a corresponding :meth:`~ase.Atoms.get_cell_lengths_and_angles`
  method.

* Galician translation of ASE's GUI.

* Two new preconditioned structure optimizers available.  See
  :mod:`ase.optimize.precon`.

* Trajectory files now contain information about the calculator and also
  information from an optimizer that wrote the trajectory.


Version 3.11.0
==============

10 May 2016: :git:`3.11.0 <../3.11.0>`.

* Special `\mathbf{k}`-points from the [Setyawan-Curtarolo]_ paper was added:
  :data:`ase.dft.kpoints.special_points`.

* New :mod:`ase.collections` module added.  Currently contains the G2 database
  of molecules and the S22 set of weakly interacting dimers and complexes.

* Moved modules:

  * ``ase.utils.eos`` moved to :mod:`ase.eos`
  * ``ase.calculators.neighborlist`` moved to :mod:`ase.neighborlist`
  * ``ase.lattice.spacegroup`` moved to :mod:`ase.spacegroup`

* The ``InfraRed`` that used to be in the ``ase.infrared`` or
  ``ase.vibrations.infrared`` modules is now called
  :class:`~ase.vibrations.Infrared` and should be imported from the
  :mod:`ase.vibrations` module.

* Deprecated modules: ``ase.structure``, ``ase.utils.geometry``,
  ``ase.utils.distance``, ``ase.lattice.surface``.  The functions from these
  modules that will create and manipulate :class:`~ase.Atoms` objects are now
  in the new :mod:`ase.build` module.  The remaining functions have been moved
  to the new :mod:`ase.geometry` module.

* The ``ase.lattice.bulk()`` function has been moved to :func:`ase.build.bulk`.

* Two new functions: :func:`~ase.geometry.cell_to_cellpar` and
  :func:`~ase.geometry.cellpar_to_cell`.

* We can now :func:`~ase.io.read` and :func:`~ase.io.write` magres files.

* :class:`~ase.neb.NEB` improvement:  calculations for molecules can now be
  told to minimize ratation and translation along the path.


Version 3.10.0
==============

17 Mar 2016: :git:`3.10.0 <../3.10.0>`.

* :ref:`old trajectory` files can no longer be used.  See :ref:`convert`.

* New iterator function :func:`ase.io.iread` for iteratively reading Atoms
  objects from a file.

* The :func:`ase.io.read` function and command-line tools can now read ``.gz``
  and ``.bz2`` compressed files.

* Two new decorators :func:`~ase.parallel.parallel_function` and
  :func:`~ase.parallel.parallel_generator` added.

* Source code moved to https://gitlab.com/ase/ase.

* Preliminary :mod:`ase.calculators.qmmm` module.

* Improved :mod:`~ase.calculators.tip3p.TIP3P` potential.

* Velocity Verlet will now work correctly with constraints.

* ASE's GUI no longer needs a special GTK-backend for matplotlib to work.
  This will make installation of ASE much simpler.

* We can now :func:`~ase.io.read` and :func:`~ase.io.write` JSV files.

* New :func:`ase.dft.kpoints.get_special_points` function.

* New :func:`ase.geometry.get_duplicate_atoms` function for finding and
  removing atoms on top of each other.

* New: A replacement :mod:`Siesta <ase.calculators.siesta>` calculator was
  implemented. It closely follows the
  :class:`ase.calculators.calculator.FileIOCalculator` class which should
  ease further development. Handling pseudopotentials, basis sets and ghost
  atoms have been made much more flexible in the new version.


Version 3.9.1
=============

21 July 2015: :git:`3.9.1 <../3.9.1>`.

* Added function for finding maximally reduced Niggli unit cell:
  :func:`ase.build.niggli_reduce`.

* Octopus interface added (experimental).


Version 3.9.0
=============

28 May 2015: :git:`3.9.0 <../3.9.0>`.

* Genetic algorithm implemented; :mod:`ase.ga`. This can be used
  for the optimization of: atomic cluster structure, materials
  properties by use of template structures. Extension to other projects
  related to atomic simulations should be straightforward.

* The ``ase.lattice.bulk`` function can now build the Wurtzite structure.

* The :class:`ase.utils.timing.Timer` was moved from GPAW to ASE.

* New :mod:`ase.db` module.

* New functions: :func:`ase.build.fcc211` and
  :func:`ase.visualize.mlab.plot`.

* New :class:`~ase.Atoms` methods:
  :meth:`ase.Atoms.get_distances()` and
  :meth:`ase.Atoms.get_all_distances()`.

* :ref:`bash completion` can now be enabled.

* Preliminary support for Python 3.

* Wrapping: new :meth:`ase.Atoms.wrap` method and
  :func:`ase.geometry.wrap_positions` function.  Also
  added ``wrap=True`` keyword argument to
  :meth:`ase.Atoms.get_scaled_positions` that can be used to turn
  off wrapping.

* New improved method for initializing NEB calculations:
  :meth:`ase.neb.NEB.interpolate`.

* New pickle-free future-proof trajectory file format added:
  :ref:`new trajectory`.

* We can now do :ref:`phase diagrams`.

* New :func:`ase.build.mx2` function for 1T and 2H metal
  dichalcogenides and friends.

* New :func:`ase.dft.bandgap.get_band_gap` function

* :class:`~ase.calculators.cp2k.CP2K` interface.


Version 3.8.0
=============

22 October 2013: :git:`3.8.0 <../3.8.0>`.

* ASE's :mod:`gui <ase.gui>` renamed from ``ag`` to ``ase-gui``.
* New :ref:`STM <stm>` module.
* Python 2.6 is now a requirement.
* The old ``ase.build.bulk`` function is now deprecated.
  Use the new one instead (:func:`ase.lattice.bulk`).
* We're now using BuildBot for continuous integration:
  https://ase-buildbot.fysik.dtu.dk/waterfall
* New interface to the JDFTx code.


Version 3.7.0
=============

13 May 2013: :git:`3.7.0 <../3.7.0>`.

* ASE's GUI can now be configured to be more friendly to visually
  impaired users: :ref:`high contrast`.

* The :class:`ase.neb.NEB` object now accepts a list of spring constants.

* *Important backwards incompatible change*: The
  :func:`ase.build.surface` function now returns a
  right-handed unit cell.

* Mopac, NWChem and Gaussian interfaces and EAM potential added.

* New :meth:`~ase.Atoms.set_initial_charges` and
  :meth:`~ase.Atoms.get_initial_charges` methods.  The
  :meth:`~ase.Atoms.get_charges` method will now ask the
  calculator to calculate the atomic charges.

* The :ref:`aep1` has been implemented and 6 ASE calculators are now
  based on the new base classes.

* ASE now runs on Windows and Mac.

* :ref:`mhtutorial` added to ASE.


Version 3.6.0
=============

24 Feb 2012: :git:`3.6.0 <../3.6.0>`.

* ASE GUI translations added, available: da_DK, en_GB, es_ES.

* New function for making surfaces with arbitrary Miller indices with
  the smallest possible surface unit cell:
  ase.build.surface()

* New ase.lattice.bulk() function.  Will replace old
  ase.build.bulk() function.  The new one will produce a more
  natural hcp lattice and it will use experimental data for crystal
  structure and lattice constants if not provided explicitly.

* New values for ase.data.covalent_radii from Cordeo *et al.*.

* New command line tool: :ref:`cli` and tests based on it:
  abinit, elk, fleur, nwchem.

* New crystal builder for ase-gui

* Van der Waals radii in ase.data

* ASE's GUI (ase-gui) now supports velocities for both graphs and coloring

* Cleaned up some name-spaces:

  * ``ase`` now contains only :class:`~ase.Atoms` and
    :class:`~ase.atom.Atom`
  * ``ase.calculators`` is now empty


Version 3.5.1
=============

24 May 2011: :git:`3.5.1 <../3.5.1>`.

* Problem with parallel vibration calculations fixed.


Version 3.5.0
=============

13 April 2011: :git:`3.5.0 <../3.5.0>`.

* Improved EMT potential:  uses a
  :class:`~ase.neighborlist.NeighborList` object and is
  now ASAP_ compatible.

* :class:`ase.optimize.BFGSLineSearch>` is now the default
  (``QuasiNewton==BFGSLineSearch``).

* There is a new interface to the LAMMPS molecular dynamics code.

* New :mod:`ase.phonons` module.

* Van der Waals corrections for DFT, see GPAW_ usage.

* New :class:`~ase.io.bundletrajectory.BundleTrajectory` added.

* Updated GUI:

  * Stability and usability improvements.
  * Povray render facility.
  * Updated expert user mode.
  * Enabled customization of colours and atomic radii.
  * Enabled user default settings via :file:`~/.ase/gui.py`.

* :mod:`Database library <ase.data>` expanded to include:

  * The s22, s26 and s22x5 sets of van der Waals bonded dimers and
    complexes by the Hobza group.
  * The DBH24 set of gas-phase reaction barrier heights by the Truhlar
    group.

* Implementation of the Dimer method.


.. _ASAP: https://wiki.fysik.dtu.dk/asap
.. _GPAW: https://gpaw.readthedocs.io/documentation/xc/vdwcorrection.html


Version 3.4.1
=============

11 August 2010: :git:`3.4.1 <../3.4.1>`.


.. toctree::

  changelog.rst
