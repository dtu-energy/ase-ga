.. module:: ase.calculators.turbomole

=========
TURBOMOLE
=========

TURBOMOLE_ is a program package for *ab initio* electronic structure calculations.
This interface integrates the TURBOMOLE code as a calculator in ASE.

.. _Turbomole: http://www.turbomole.com/


Setting up the environment
==========================

The TURBOMOLE package must be installed to use it with ASE. All modules and
scripts from the TURBOMOLE packages must be available in $PATH and the variable
$TURBODIR must be set. More information on how to install TURBOMOLE and to set
up the environment can be found in the manual or the tutorial at
the `web site`_.

.. _web site: http://www.turbomole-gmbh.com/turbomole-manuals.html

Using the calculator
====================

Python interface
----------------

The constructor method has only keyword arguments that can be specified in any
order. The list of accepted parameters with their types and default values is
provided in the section "Parameters" below.

The following example demonstrates how to construct a Turbomole calculator
object for a single-point energy calculation of a neutral singlet
system:

.. code:: python

  from ase.calculators.turbomole import Turbomole
  calc = Turbomole(multiplicity=1)

The selection of the method will be according to the default parameter values
(see below), i.e. in this case DFT with b-p functional and the def-SV(P) basis
set. After this the calculator can be associated with an existing Atoms object

.. code:: python

  atoms.calc = calc

The recommended methods to access parameters and properties are the getter
methods, i.e. these ones starting with *get*. The calculations then are
triggered according to the principle of lazy evaluation, i.g.:

.. code:: python

  energy = atoms.get_potential_energy()
  print(energy)

Alternatively all calculations necessary to perform a task (see ``task``
parameter below) can be explicitly started with the ``calculate()`` method:

.. code:: python

  calc.calculate(atoms)

The getter methods (see below) check for convergence and eventually return
``None`` or an exception if the calculation has not converged. If the
properties are read using the Turbomole object attributes then the convergence
must be checked with:

.. code:: python

  assert calc.converged

If the user wishes to use the input files (such as the control file) generated
by module ``define`` before (or without) an actual calculation starts, the
``initialize()`` method has to be called explicitly after constructing the
calculator and associating it with an atoms object, e.g.:

.. code:: python

    from ase.build import molecule
    from ase.calculators.turbomole import Turbomole
    mol = molecule('C60')
    params = {
        'use resolution of identity': True,
        'total charge': -1,
        'multiplicity': 2
    }
    calc = Turbomole(**params)
    mol.calc = calc
    calc.initialize()




Command-line interface
----------------------

The command-line interface has limited capability. For example the keyword
``task`` is not effective due to the specific way the methods are called by
``ase-run``. This example shows how to run a single-point DFT calculation of
water with the PBE functional and with geometry taken from the database::

  ase-build H2O | ase-run turbomole --parameters="multiplicity=1,density functional=pbe"

Using the calculation output a second geometry optimization calculation with the
BFGS optimizer from ASE can be started using the ``restart`` keyword::

  ase-build H2O | ase-run turbomole --parameters="restart=True" -f 0.02


Reading output
==============

Properties
----------

The implemented properties are described in the following table.

================== ======== ======================= =========== ==================
**Property**       **Type** **Getter method**       **Storage** **Task**
================== ======== ======================= =========== ==================
total energy       float    get_potential_energy(), e_total     any task
                            get_property('energy')
forces             np.array get_forces(),           forces      gradient
                            get_property('forces')
dipole moment      np.array get_dipole_moment(),    dipole      any task
                            get_property('magmom')
charges            np.array get_charges(),
                            get_property('charges') charges     any task
<S\ :sup:`2`\ >    float    get_results             results     any task
normal modes       list     get_results             results     frequencies
mode frequencies   list     get_results             results     frequencies
gradient           list     get_results             results     gradient, optimize
hessian            list     get_results             results     frequencies
molecular orbitals list     get_results             results     any task
occupancies        list     get_results             results     any task
================== ======== ======================= =========== ==================

Metadata
--------

Additionally, some useful information can be read with the calculator using the
functions ``read_version()``, ``read_datetime()``, ``read_runtime()``,
``read_hostname()``. Then the respective data can be retrieved using the
*version*, *datetime*, *runtime* and *hostname* attributes. Example:

.. code:: python

  calc.read_runtime()
  print(calc.runtime)


Restart mode
------------

The restart mode can be used either to start a calculation from the data left
from previous calculations or to analyze or post-process these data. The
previous run may have been performed without ASE but the working directory of
the job should contain the control file and all files referenced in it. In
addition, the standard output will be searched in files beginning with *job.*
and ending with *.out* but this is optional input, mainly to extract job
datetime, runtimes, hostname and TURBOMOLE version. After constructing the
calculator object (where *params* dictionary is optional):

.. code:: python

  calc = Turbomole(restart=True, **params)

the data left from the previous calculations can be queried, for example:

.. code:: python

  from ase.visualize import view
  view(calc.atoms)
  print(calc.converged)
  print(calc.get_potential_energy())

A previous calculation may have crashed or not converged. Also in these cases
all data that is available will be retrieved but the ``calc.converged`` will
be set to ``False``. The calculation can be continued without any parameter
modifications (for example if it has exceeded the job maximum run time and was
interrupted) or with better convergence parameters specified in ``params``
dictionary. Finally, another calculation task can be started beginning
from the data left from a converged previous one, specifying a new ``task``
parameter:

.. code:: python

  calc = Turbomole(restart=True, task='gradient', **params)


Caveat about using the restart mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When restart mode is set, the calculator reads and parses all relevant data
groups from the *control* file. The format of the `control` file, especially
of some data groups, does not allow to guarantee that the parsing is always
successful and correct. The parsing process may fail and raise an exception
or the parameters read might be wrong. If you encounter errors due to restart
mode please file an issue.


Policies for files in the working directory
-------------------------------------------

* When the calculator is constructed in restart mode (i.e. ``restart=True``)
  and with no other parameters, then no files will be created, deleted or
  modified in the working directory.

* When the calculator is created in normal (i.e. ``restart=False``) mode then
  all TURBOMOLE related files found in the working directory will be deleted.

* When the calculator is created with ``restart=True`` and other parameters,
  the *control* file might be modified. In particular, if ``define_str``,
  ``control_input`` or ``control_kdg`` are specified or ``initialize()``
  is called then the *control* file will be modified.

* When ``calculate()``, ``get_potential_energy()``, ``get_forces()`` etc. are
  called in restart mode, the *control* file will be modified if the previous
  calculation has not converged.

* When an *atoms* object is associated with the calculator or any calculator
  method is called with an *atoms* object specified, then the calculator will
  be reset and all TURBOMOLE related files found in the working directory will
  be deleted if *atoms* is different (tol=1e-2) from the internal *atoms* object or
  if internal coordinates are used and the internal and the supplied *atoms*
  positions are different (tol=1e-13). The *coord* file will be changed only
  if the *atoms* positions are different (tol=1e-13).


Parameters
==========

The following table provides a summary of all parameters and their default
values.

================================ ======== =========== ============= ==============
**Name**                         **Type** **Default** **Units**     **Updateable**
================================ ======== =========== ============= ==============
                         restart  bool    False       None            True
                      define_str   str    None        None            True
                     control_kdg  list    None        None            True
                   control_input  list    None        None            True
                 reset_tolerance float    1e-2        Angstrom        True
         automatic orbital shift float          0.1             eV          True
                  basis set name   str    def-SV(P)           None         False
      closed-shell orbital shift float         None             eV          True
         damping adjustment step float         None           None          True
             density convergence float         None           None          True
              density functional   str          b-p           None          True
              energy convergence float         None             eV          True
                         esp fit   str         None           None          True
          fermi annealing factor float         0.95           None          True
         fermi final temperature float          300         Kelvin          True
   fermi homo-lumo gap criterion float          0.1             eV          True
       fermi initial temperature float          300         Kelvin          True
        fermi stopping criterion float        0.001             eV          True
               force convergence float         None    eV/Angstrom          True
geometry optimization iterations   int         None           None          True
                       grid size   str           m3           None          True
                    ground state  bool         True           None         False
                 initial damping float         None           None          True
                   initial guess  None          eht           None         False
                 minimal damping float         None           None          True
                    multiplicity   int         None           None         False
     non-automatic orbital shift  bool        False           None          True
               numerical hessian  dict         None           None          True
                     point group   str           c1           None         False
                       ri memory   int         1000       Megabyte          True
          scf energy convergence float         None             eV          True
                  scf iterations   int           60           None          True
                            task   str       energy           None          True
                           title   str           ''           None         False
                    total charge   int            0           None         False
               transition vector   int         None           None          True
                             uhf  bool         None           None         False
           use basis set library  bool         True           None         False
                         use dft  bool         True           None         False
              use fermi smearing  bool        False           None          True
         use redundant internals  bool        False           None         False
      use resolution of identity  bool        False           None         False
================================ ======== =========== ============= ==============

The attribute ``Updateable`` specifies whether it is possible to change a
parameter upon restart. The ``restart`` keyword tells the calculator whether to
restart from a previous calculation. The optional ``define_str`` is a string of
characters that would be entered in an interactive session with module ``define``,
i.e. this is the stdin for running module ``define``. The ``control_kdg`` is an
optional list of data groups in control file to be deleted after running module
``define`` and ``control_input`` is an optional list of data groups to be added
to control file after running module ``define``.

If the Atoms object is updated via ``set_atoms()`` method, a check for the changes
is performed and if the changes in positions are larger than a tolerance
``reset_tolerance`` then the calculator is reset, the working directory is purged
and module ``define`` is called. In order to control this behavior the user may
choose a custom value for ``reset_tolerance``.

The parameter ``initial guess`` can be either the strings *eht* (extended
Hückel theory) or *hcore* (one-electron core Hamiltonian) or a dictionary
*{'use': '<path/to/control>'}* specifying a path to a control file with the
molecular orbitals that should be used as initial guess.

If ``numerical hessian`` is defined then the force constant matrix will be
computed numerically using the script NumForce. The keys can be *'central'*
indicating use of central differences (type *bool*) and *'delta'* specifying
the coordinate displacements in Angstrom (type *float*).

While ``task`` can be set to ``"optimize"`` to perform a geometry optimization
using Turbomole's own relaxation algorithms, doing so directly is discouraged.
Instead, the calculator's ``get_optimizer()`` method should be called to obtain
a ``TurbomoleOptimizer`` which can be used like any other ASE
:mod:`Optimizer <ase.optimize>`. An :ref:`example <turbomole_optimizer_example>`
is given below.

Some parameter names contain spaces. This means that the preferred way to pass
the parameters is to construct a dictionary, for example:

.. code:: python

  params = {'use resolution of identity': True,
            'ri memory': 2000,
            'scf iterations': 80,
            'force convergence': 0.05}
  calc = Turbomole(**params)

Using the ``todict()`` method, the parameters of an existing Turbomole calculator
object can be stored in a flat dictionary and then re-used to create a
new Turbomole calculator object:

.. code:: python

  params = calc.todict()
  new_calc = Turbomole(**params)

This is especially useful if the *calc* object has been created in restart
mode or retrieved from a database.


Examples
========

Single-point energy calculation
-------------------------------

This script calculates the total energy of H2:

:git:`ase/test/calculator/turbomole/test_turbomole_H2.py`.

Nudged elastic band calculation
-------------------------------

The example demonstrates a proton transfer barrier calculation in H3O2-:

:git:`ase/test/calculator/turbomole/test_turbomole_h3o2m.py`.

Single-point gradient calculation of Au13-
------------------------------------------

This script demonstrates the use of the restart option.

:git:`ase/test/calculator/turbomole/test_turbomole_au13.py`.


.. _turbomole_optimizer_example:

Geometry optimization using TurbomoleOptimizer (recommended)
------------------------------------------------------------

:git:`ase/test/calculator/turbomole/test_turbomole_optimizer.py`.

Geometry optimization and normal mode analysis for H2O
------------------------------------------------------

:git:`ase/test/calculator/turbomole/test_turbomole_h2o.py`.


.. _turbomole qmmm:

QMMM simulation
---------------

The following example demonstrates how to use the Turbomole calculator in simple
and explicit QMMM simulations on the examples of a water dimer partitioned into
an MM and a QM region.

:git:`ase/test/calculator/turbomole/test_turbomole_qmmm.py`.

The MM region is treated within a TIP3P model in the MM calculator and as an
array of point charges in the QM calculation. The interaction between the QM
and MM regions, used in the explicit QMMM calculator, is of Lennard-Jones type.

The point charge embedding functionality of the Turbomole calculator can also be
used without QMMM calculators if the ``embed()`` method is called with a
specification of the point charges and their positions in which to embed the
QM system:

.. code:: python

    from ase.collections import s22
    from ase.calculators.turbomole import Turbomole

    params = {'esp fit': 'kollman', 'multiplicity': 1}
    dimer = s22['Water_dimer']
    qm_mol = dimer[0:3]
    calc = Turbomole(**params)
    qm_mol.calc = calc
    calc.embed(
        charges=[-0.76, 0.38,  0.38],
        positions=dimer.positions[3:6]
    )
    print(qm_mol.get_potential_energy())
    print(qm_mol.get_forces())
    print(qm_mol.get_charges())

A more elaborated version of the latter example is used in the test script:

:git:`ase/test/calculator/turbomole/test_turbomole_2h2o.py`.


Deprecated, non-implemented and unsupported features
====================================================

Deprecated but still accepted parameters
----------------------------------------

==================== ======== ======================== =========================
Name                 Type     Default value            Description
==================== ======== ======================== =========================
``calculate_energy`` ``str``  ``dscf``                 module name for energy
                                                       calculation
``calculate_forces`` ``str``  ``grad``                 module name for forces
                                                       calculation
``post_HF``          ``bool``  ``False``               post Hartree-Fock format
                                                       for energy reader
==================== ======== ======================== =========================


Not implemented parameters
--------------------------

The following table includes parameters that are planned but not implemented yet.

================================ ======= ========== =============== ==========
Name                             Type    Default    Units           Updateable
================================ ======= ========== =============== ==========
            basis set definition  dict        None          None         False
                   excited state  bool       False          None         False
                           label   str        None          None         False
        number of excited states   int        None          None         False
         optimized excited state   int        None          None         False
                            rohf  bool        None          None         False
================================ ======= ========== =============== ==========


Unsupported methods and features
--------------------------------

The following methods and features are supported in TURBOMOLE but currently not
in the ASE Turbomole calculator:

* MP2 and coupled-cluster methods (modules mpgrad, rimp2, ricc2)
* Excited state calculations (modules escf, egrad)
* Molecular dynamics (modules mdprep, uff)
* Solvent effects (COSMO model)
* Global optimization (module haga)
* Property modules (modules freeh, moloch)
* Point groups other than C1 (see not implemented parameters)
* Restricted open-shell Hartree-Fock (see not implemented parameters)
* Per-element and per-atom basis set specifications (see not implemented parameters)
* Explicit basis set specification (see not implemented parameters)
