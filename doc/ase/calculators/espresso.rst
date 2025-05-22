.. module:: ase.calculators.espresso

========
Espresso
========

.. image:: ../../static/espresso.png

`Quantum ESPRESSO <http://www.quantum-espresso.org>`_ (QE) is an integrated
suite of Open-Source computer codes for electronic-structure calculations and
materials modeling at the nanoscale. It is based on density-functional
theory, plane waves, and pseudopotentials.
The ASE calculator is an interface to the ``pw.x`` executable, however,
input/output operations can be managed for other executables (see the Input/Output section below). Users interested to run calculation with other executables,
can have a look in the ase ecosystem pages for other packages.

Calculator
==========

Basic usage
-----------

Any calculation will need pseudopotentials for the elements involved. The
default directory for the pseudopotential files can be set with the ``pseudo_dir``
parameter in the ASE configuration.  Use ``ase info --calculators espresso``
Pseudopotential files can be passed
to the calculator as a dictionary, where the keys are the element and values
are the pseudopotential file names. The pseudopotential file names must be
in the directory specified by the ``pseudo_dir`` parameter.

.. code-block:: python

    pseudopotentials = {"Na": "na_pbe_v1.5.uspp.F.UPF", "Cl": "cl_pbe_v1.4.uspp.F.UPF"}

A recommended list of pseudopotentials is the SSSP curated set, which can be
found on `Materials Cloud <https://www.materialscloud.org
/discover/sssp/table/efficiency>`_.
With the all of the above in mind, a simple calculation can be set up like so:

.. code-block:: python

  from ase.build import bulk
  from ase.calculators.espresso import Espresso, EspressoProfile
  from ase.optimize import LBFGS

  rocksalt = bulk('NaCl', crystalstructure='rocksalt', a=6.0)

  # Pseudopotentials from SSSP Efficiency v1.3.0
  pseudopotentials = {'Na': 'na_pbe_v1.5.uspp.F.UPF', 'Cl': 'cl_pbe_v1.4.uspp.F.UPF'}

  # Optionally create profile to override paths in ASE configuration:
  profile = EspressoProfile(
      command='/path/to/pw.x', pseudo_dir='/path/to/pseudopotentials'
  )

  calc = Espresso(profile=profile, pseudopotentials=pseudopotentials)

  rocksalt.calc = calc

  rocksalt.get_potential_energy()  # This will run a single point calculation

  opt = LBFGS(rocksalt)

  opt.run(fmax=0.005)  # This will run a geometry optimization using ASE's LBFGS algorithm

  # Print lattice constant...
  print((8 * rocksalt.get_volume() / len(rocksalt)) ** (1.0 / 3.0))


The above example will run a single point calculation and a geometry optimization
using default QE settings. In the next sections, we will describe how to
customize the QE settings.

Parameters
----------

The calculator will interpret most documented options for ``pw.x``
found at the `PW input description <https://www.quantum-espresso.org/
Doc/INPUT_PW.html>`_.

Parameters must be passed in the ``input_data`` parameter, which is a dictionary
, previously accepted kwargs keywords are now deprecated. The ``input_data`` dictionary
can be used both in flat or nested format. In the flat format, keywords
will be put into the correct section of the input file. In the nested format, the
``input_data`` dictionary will be used as is. Currently the nested format can be used
to specify exotic input sections manually such as ``&FCP`` or ``&RISM``.

All parameters must be given in QE units, usually Ry or atomic units
in line with the documentation. ASE does not add any defaults over the
defaults of QE.

.. code-block:: python

  input_data = {
      'system': {'ecutwfc': 60, 'ecutrho': 480},
      'disk_io': 'low',  # Automatically put into the 'control' section
  }

  calc = Espresso(
      profile=profile
      pseudopotentials=pseudopotentials,
      tstress=True,  # deprecated, put in input_data
      tprnfor=True,  # deprecated, put in input_data
      input_data=input_data,
  )

Any ``FixAtoms`` or ``FixCartesian`` constraints are converted to Espresso
constraints. Some parameters are used by ASE and have additional meaning:

 * ``kpts``, is used to specify the k-point sampling
   * If ``kpts`` is a tuple (or list) of 3 integers ``kpts=(int, int, int)``, it is interpreted  as the dimensions of a Monkhorst-Pack grid.
   * If ``kpts`` is set to ``None``, only the Γ-point will be included and QE will use routines optimized for Γ-point-only calculations. Compared to Γ-point-only calculations without this optimization (i.e. with ``kpts=(1, 1, 1)``), the memory and CPU requirements are typically reduced by half.
   * If ``kpts`` is a dict, it will either be interpreted as a path in the Brillouin zone (see ase.dft.kpoints.bandpath) if it contains the 'path' keyword, otherwise it is converted to a Monkhorst-Pack grid (see ase.calculators.calculator.kpts2sizeandoffsets)
 * ``koffset=(0, 0, 0)`` set to 0 or 1 to displace the kpoint grid by a half
   cell in that direction. Also accepts ``True`` and ``False``.
 * ``kspacing=0.1`` sets the minimum distance between kpoints in reciprocal
   space.
 * ``nspin`` does not typically need to be specified by the user. If any atom has a magnetic moment, it is turned on automatically.

Parallelism
-----------

The above ``calc`` object will run the ``pw.x`` executable using the command
specified by the ASE configuration.  Parallelism
can be customized by manually creating an ``EspressoProfile``.
Currently there is no way to specify the parallelization keywords for QE, such
as ``-nk``. This is a current limitation of the ASE calculator interface.

.. code-block:: python

  profile = EspressoProfile(
      command='mpiexec -n 16 /path/to/pw.x',
      pseudo_dir='/path/to/pseudopotentials',
  )

  # This will run 'mpirun -np 16 -v /path/to/pw.x'

Input/Output
============

ASE can read and write input files for the espresso executables. The procedure
depends on the executable. Each subsection will describe the procedure for each
executable.

pw.x
----

The ``pw.x`` executable is the main executable of the Quantum Espresso suite. ASE
can fully read and write input files for this executable. This must be done by
using the main :func:`ase.io.write` and :func:`ase.io.read` functions. The

These two main functions are wrappers around the :func:`ase.io.espresso.write_espresso_in` and :func:`ase.io.espresso.read_espresso_out` functions. The former will write an input file for the ``pw.x`` executable, while the latter will read an input file and return Atoms object(s) with available results.

.. autofunction:: ase.io.espresso.write_espresso_in
  :noindex:

**Example:**

.. code-block:: python

  from ase.io.espresso import write_espresso_in

  input_data = {
      'calculation': 'relax',
      'restart_mode': 'from_scratch',
      'tprnfor': True,
      'etot_conv_thr': 1e-5,
      'forc_conv_thr': 1e-4,
      'ecutwfc': 60,
      'ecutrho': 480,
      'input_dft': 'rpbe',
      'vdw_corr': 'dft-d3',
      'occupations': 'smearing',
      'degauss': 0.01,
      'smearing': 'cold',
      'conv_thr': 1e-8,
      'mixing_mode': 'local-TF',
      'mixing_beta': 0.35,
      'diagonalization': 'david',
      'ion_dynamics': 'bfgs',
      'bfgs_ndim': 6,
      'startingwfc': 'random',
  }  # This flat dictionary will be converted to a nested dictionary where, for example, "calculation" will be put into the "control" section

  write('pw.in', atoms, input_data=input_data, format='espresso-in')

  # Hubbard is not implemented in the write_espresso_in function, we can add it manually
  additional_cards = ['HUBBARD (ortho-atomic)', 'U Mn-3d 5.0', 'U Ni-3d 6.0']

  write(
      'pw_hubbard.in',
      atoms,
      input_data=input_data,
      additional_cards=additional_cards,
      format='espresso-in',
  )

.. autofunction:: ase.io.espresso.read_espresso_out
  :noindex:

ph.x
----

As of January 2024, ``ph.x`` has custom read and write functions in ASE. The function :meth:`~ase.io.espresso.write_espresso_ph` can be used to write a Fortran namelist file for ``ph.x``.

.. autofunction:: ase.io.espresso.write_espresso_ph

**Example:**

.. code-block:: python

  from ase.io.espresso import write_espresso_ph, read_espresso_ph

  input_data = {
    'tr2_ph': 1.0e-12,
    'prefix': 'pwscf',
    'verbosity': 'high',
    'ldisp': True,
    'qplot': True,
    'alpha_mix(1)': 0.1,
  }

  qpts = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5), (0.5, 0.5, 0.0), (0.0, 0.0, 0.5)]

  write_espresso_ph('input_ph.in', input_data, qpts=qpts) # Will automatically be built to nested format

After running the calculation, the output can be read with the :meth:`~ase.io.espresso.read_espresso_ph` method.

.. autofunction:: ase.io.espresso.read_espresso_ph

**Example:**

.. code-block:: python

  results = read_espresso_ph('ph.out')

  print(results[0]['qpoints']) # Will probably be (0, 0, 0)
  print(results[0]['freq']) # Eigenvalues of the dynamical matrix at the qpoint

Other binaries
--------------

Since January 2024, it is possible to use ASE to write input files for other
QE executables.

The :func:`ase.io.espresso.write_fortran_namelist` method
will write a fortran namelist file, which can be read by QE executables.
The list of currently implemented executable is available in
:mod:`ase.io.espresso_namelist.keys`

.. autofunction:: ase.io.espresso.write_fortran_namelist

**Example:**

.. code-block:: python

  from ase.io.espresso import Namelist
  from ase.io.espresso import write_fortran_namelist

  # input_data for pp.x, built automatically
  input_data = {'outdir': '/path/to/output', 'prefix': 'prefix', 'verbosity': 'high'}

  input_data = Namelist(input_data)
  input_data.to_nested(binary='pp')

  write_fortran_namelist('input_pp.in', input_data)

  # Alternatively, the input_data can be written manually
  # without the need to specify a binary

  input_data = {
    'environ': {
      'environ_type': 'water',
    ...} # In this case, the ``input_data`` is a nested dictionary
    # since ASE cannot guess the input sections for exotic executables

  additional_cards = [
        'EXTERNAL_CHARGES (bohr)',
        '-0.5 0. 0. 25.697 1.0 2 3',
        '-0.5 0. 0. 20.697 1.0 2 3'
  ] # This will be added at the end of the file

  write_fortran_namelist('environ.in', input_data)

The :meth:`~ase.io.espresso.read_fortran_namelist` method will read a fortran namelist file and
return a tuple, where the first element is the ``input_data`` dictionary and the
second element is a list of additional cards at the end of the file.

.. autofunction:: ase.io.espresso.read_fortran_namelist

**Example:**

.. code-block:: python

  from ase.io.espresso import read_fortran_namelist

  input_data, additional_cards = read_fortran_namelist('input_pp.in')

  print(input_data)
  print(additional_cards)

You can use both of these methods to jungle between your favorite QE executables.

Namelist Class
==============

The Namelist class is a simple class to handle fortran namelist files. It derived from the dict class and is case case-insensitive.

.. autoclass:: ase.io.espresso.Namelist
      

Espresso Calculator Class
=========================

.. autoclass:: ase.calculators.espresso.Espresso

