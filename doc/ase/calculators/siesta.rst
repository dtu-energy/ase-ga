.. module:: ase.calculators.siesta

======
SIESTA
======

Introduction
============

SIESTA_ is a density-functional theory code for very large systems
based on atomic orbital (LCAO) basis sets.


.. _SIESTA: https://siesta-project.org/siesta/



Environment variables
=====================

The environment variable :envvar:`ASE_SIESTA_COMMAND` must hold the command
to invoke the siesta calculation. The variable must be a string 
where ``PREFIX.fdf``/``PREFIX.out`` are the placeholders for the 
input/output files. This variable allows you to specify serial or parallel 
execution of SIESTA.
Examples: ``siesta < PREFIX.fdf > PREFIX.out`` and
``mpirun -np 4 /bin/siesta-5.0 < PREFIX.fdf > PREFIX.out``.

A default directory holding pseudopotential files :file:`.psml/.psf` can be
defined to avoid defining this every time the calculator is used.
This directory can be set by the environment variable
:envvar:`SIESTA_PP_PATH`.

Set both environment variables in your shell configuration file:

.. highlight:: bash

::

  $ export ASE_SIESTA_COMMAND="siesta < PREFIX.fdf > PREFIX.out"
  $ export SIESTA_PP_PATH=$HOME/mypps

.. highlight:: python

Alternatively, the path to the pseudopotentials can be given in
the calculator initialization. Listed below all parameters 
related to pseudopotential control.

===================== ========= ============= =====================================
keyword               type      default value description
===================== ========= ============= =====================================
``pseudo_path``       ``str``   ``None``      Directory for pseudopotentials to use
                                              None means using $SIESTA_PP_PATH
``pseudo_qualifier``  ``str``   ``None``      String for picking out specific type
                                              type of pseudopotentials. Giving
                                              ``example`` means that
                                              ``H.example.psml`` or
                                              ``H.example.psf`` will be used. None
                                              means that the XC.functional keyword
                                              is used, e.g. ``H.lda.psf``
``symlink_pseudos``   ``bool``  ``---``       Whether pseudos will be sym-linked 
                                              into the execution directory. If 
                                              False they will be copied instead.
                                              Default is True on Unix and False on
                                              Windows.
===================== ========= ============= =====================================


SIESTA Calculator
=================

These parameters are set explicitly and overrides the native values if different.

================ ========= =================== =====================================
keyword          type      default value       description
================ ========= =================== =====================================
``label``        ``str``   ``'siesta'``        Name of the output file
``mesh_cutoff``  ``float`` ``200*Ry``          Mesh cut-off energy in eV
``xc``           ``str``   ``'LDA'``           Exchange-correlation functional.
                                               Corresponds to either XC.functional
                                               or XC.authors keyword in SIESTA
``energy_shift`` ``float`` ``100 meV``         Energy shift for determining cutoff
                                               radii
``kpts``         ``list``  ``[1,1,1]``         Monkhorst-Pack k-point sampling
``basis_set``    ``str``   ``DZP``             Type of basis set ('SZ', 'DZ', 'SZP',
                                               'DZP')
``spin``         ``str``   ``'non-polarized'`` The spin approximation used, must be
                                               either ``'non-polarized'``, 
                                               ``'collinear'``, ``'non-collinear'``
                                               or ``'spin-orbit'``.
``species``      ``list``  ``[]``              A method for specifying the basis set  
                                               for some atoms.
================ ========= =================== =====================================

Most other parameters are set to the default values of the native interface.

Extra FDF parameters
====================

The SIESTA code reads the input parameters for any calculation from a
:file:`.fdf` file. This means that you can set parameters by manually setting
entries in this input :file:`.fdf` file. This is done by the argument:

>>> Siesta(fdf_arguments={'variable_name': value, 'other_name': other_value})

For example, the ``DM.MixingWeight`` can be set using

>>> Siesta(fdf_arguments={'DM.MixingWeight': 0.01})

The explicit fdf arguments will always override those given by other
keywords, even if it breaks calculator functionality.
The complete list of the FDF entries can be found in the official `SIESTA
manual`_.

.. _SIESTA manual:  https://siesta-project.org/siesta/Documentation/

Example
=======

Here is an example of how to calculate the total energy for bulk Silicon,
using a double-zeta basis generated by specifying a given energy-shift:

>>> from ase import Atoms
>>> from ase.calculators.siesta import Siesta
>>> from ase.units import Ry
>>>
>>> a0 = 5.43
>>> bulk = Atoms('Si2', [(0, 0, 0),
...                      (0.25, 0.25, 0.25)],
...              pbc=True)
>>> b = a0 / 2
>>> bulk.set_cell([(0, b, b),
...                (b, 0, b),
...                (b, b, 0)], scale_atoms=True)
>>>
>>> calc = Siesta(label='Si',
...               xc='PBE',
...               mesh_cutoff=200 * Ry,
...               energy_shift=0.01 * Ry,
...               basis_set='DZ',
...               kpts=[10, 10, 10],
...               fdf_arguments={'DM.MixingWeight': 0.1,
...                              'MaxSCFIterations': 100},
...               )
>>> bulk.calc = calc
>>> e = bulk.get_potential_energy()

Here, the only input information on the basis set is, that it should
be double-zeta (``basis='DZP'``) and that the confinement potential
should result in an energy shift of 0.01 Rydberg (the
``energy_shift=0.01 * Ry`` keyword). Sometimes it can be necessary to specify
more information on the basis set.

Defining Custom Species
=======================
Standard basis sets can be set by the keyword ``basis_set`` directly, but for
anything more complex than one standard basis size for all elements,
a list of ``species`` must be defined. Each specie is identified by atomic
element and the tag set on the atom.

For instance if we wish to investigate a H2 molecule and put a ghost atom
(the basis set corresponding to an atom but without the actual atom) in the middle
with a special type of basis set you would write:

>>> from ase.calculators.siesta.parameters import Specie, PAOBasisBlock
>>> from ase import Atoms
>>> from ase.calculators.siesta import Siesta
>>> atoms = Atoms(
...     '3H',
...     [(0.0, 0.0, 0.0),
...      (0.0, 0.0, 0.5),
...      (0.0, 0.0, 1.0)],
...     cell=[10, 10, 10])
>>> atoms.set_tags([0, 1, 0])
>>>
>>> basis_set = PAOBasisBlock(
... """1
... 0  2 S 0.2
... 0.0 0.0""")
>>>
>>> siesta = Siesta(
...     species=[
...         Specie(symbol='H', tag=None, basis_set='SZ'),
...         Specie(symbol='H', tag=1, basis_set=basis_set, ghost=True)])
>>>
>>> atoms.calc = siesta

When more species are defined, species defined with a tag has the highest priority.
General species with ``tag=None`` has a lower priority.
Finally, if no species apply
to an atom, the general calculator keywords are used.


Pseudopotentials
================

Siesta requires for pseudopotential files in the ``.psml`` or ``.psf`` formats. 
Up-to-date information on Siesta compatible pseudopotential formats and databases 
can be accessed on the `pseudopotentials section`_ of the siesta-project website.

* Curated pseudopotential databases in ``.psml`` format are available at PseudoDojo_. These are supported by SIESTA 5.0 and later versions.
* Optimized GGA–PBE pseudos in ``.psf`` format and DZP basis sets for some common elements are also available from the `SIMUNE`_ website.
* PSML pseudopotentials can be generated with the `ONCVPSP generator`_.

.. _pseudopotentials section: https://siesta-project.org/siesta/Documentation/Pseudopotentials/

.. _PseudoDojo: http://www.pseudo-dojo.org/

.. _SIMUNE:  https://www.simuneatomistics.com/siesta-toolkit/siesta-pseudos-and-basis-database/

.. _ONCVPSP generator: https://github.com/oncvpsp/oncvpsp


Species can also be used to specify pseudopotentials:

>>> specie = Specie(symbol='H', tag=1, pseudopotential='H.example.psml')

When specifying the pseudopotential in this manner, both absolute
and relative paths can be given.
Relative paths are interpreted as relative to the set 
pseudopotential path.

Restarting from an old Calculation
==================================

If you want to rerun an old SIESTA calculation, whether made using the ASE
interface or not, you can set the keyword ``restart`` to the siesta ``.XV``
file. The keyword ``ignore_bad_restart`` (True/False) will decide whether
a broken file will result in an error(False) or the whether the calculator
will simply continue without the restart file.

Choosing the coordinate format
==============================
If you are mainly using ASE to generate SIESTA files for relaxation with native
SIESTA relaxation, you may want to write the coordinates in the Z-matrix format
which will best allow you to use the advanced constraints present in SIESTA.

======================= ========= ============= =====================================
keyword                 type      default value description
======================= ========= ============= =====================================
``atomic_coord_format`` ``str``   ``'xyz'``     Choose between ``'xyz'`` and 
                                                ``'zmatrix'`` for the format that 
                                                coordinates will be written in.
======================= ========= ============= =====================================

Siesta Calculator Class
=======================

.. autoclass:: ase.calculators.siesta.siesta.Siesta


Excited states calculations
===========================

The `PyNAO <https://mbarbrywebsite.ddns.net/pynao/doc/html/>`_ code can be used
to access excited state properties after having obtained the ground state
properties with SIESTA. PyNAO allows to perform

* Time Dependent Density Functional Theory (TDDFT) calculations
* GW approximation calculations
* Bethe-Salpeter equation (BSE) calculations.

Example of code to calculate polarizability of CH4 molecule::

  from ase.calculators.siesta.siesta_lrtddft import SiestaLRTDDFT
  from ase.build import molecule
  import numpy as np

  # Define the systems
  ch4 = molecule('CH4')
  lrtddft = SiestaLRTDDFT(label="siesta", xc_code='LDA,PZ')

  # run DFT with siesta
  lrtddft.get_ground_state(ch4)

  # Run TDDFT with PyNAO
  freq = np.arange(0.0, 25.0, 0.5)
  pmat = lrtddft.get_polarizability(freq)

  import matplotlib.pyplot as plt

  plt.plot(freq, pmat[0, 0, :].imag)
  plt.show()

Raman Calculations with SIESTA and PyNAO
========================================

It is possible to calculate the Raman spectra with SIESTA and PyNAO using the
Raman function of the vibration module::

  from ase.calculators.siesta.siesta_lrtddft import RamanCalculatorInterface
  from ase.calculators.siesta import Siesta
  from ase.vibrations.raman import StaticRamanCalculator
  from ase.vibrations.placzek import PlaczekStatic
  from ase.build import molecule

  n2 = molecule('N2')

  # enter siesta input
  n2.calc = Siesta(
      basis_set='DZP',
      fdf_arguments={
          'COOP.Write': True,
          'WriteDenchar': True,
          'XML.Write': True})

  name = 'n2'
  pynao_args = dict(label="siesta", jcutoff=7, iter_broadening=0.15,
                    xc_code='LDA,PZ', tol_loc=1e-6, tol_biloc=1e-7)
  rm = StaticRamanCalculator(n2, RamanCalculatorInterface, name=name,
                             delta=0.011, exkwargs=pynao_args)
  rm.run()

  Pz = PlaczekStatic(n2, name=name)
  e_vib = Pz.get_energies()
  Pz.summary()

Further Examples
================
See also ``ase/test/calculators/siesta/lrtddft`` for further examples
on how the calculator can be used.

Siesta lrtddft Class
====================

.. autoclass:: ase.calculators.siesta.siesta_lrtddft.SiestaLRTDDFT

Siesta RamanCalculatorInterface Calculator Class
================================================

.. autoclass:: ase.calculators.siesta.siesta_lrtddft.RamanCalculatorInterface
