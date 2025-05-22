.. module:: ase.calculators.nwchem

======
NWChem
======

`NWChem <https://nwchemgit.github.io>`_ is a computational chemistry code
based on gaussian basis functions or plane-waves.


Setup
=====

.. highlight:: bash

You first need to install a working copy of NWChem for ASE to call;
follow the instructions on the `NWChem website <https://nwchemgit.github.io>`_.

The default command that ASE will use to start NWChem is
``nwchem PREFIX.nwi > PREFIX.nwo``. You can change this command by setting the
environment variable :envvar:`ASE_NWCHEM_COMMAND`. (For example, add a line
to your ``.bashrc`` with ``export ASE_NWCHEM_COMMAND="my new command"``.)

The default command will only allow you to run NWChem on a single core. To
run on multiple processors you will need to specify your MPI (or similar)
command, which typically requires telling the MPI command the number of tasks
to dedicate to the process. An example command to allow multiprocessing is
``mpirun -n $SLURM_NTASKS nwchem PREFIX.nwi > PREFIX.nwo``, for the SLURM
queueing system. If you use a different queueing system replace
``$SLURM_NTASKS`` with the appropriate variable, such as ``$PBS_NP``.


Examples
========

Here is a command line example of how to optimize the geometry of a
water molecule using the PBE density functional::

    $ ase build H2O | ase run nwchem -p xc=PBE -f 0.02
    Running: H2O
    LBFGS:   0  09:58:54    -2064.914841       1.9673
    LBFGS:   1  09:58:55    -2064.976691       0.1723
    LBFGS:   2  09:58:55    -2064.977120       0.0642
    LBFGS:   3  09:58:55    -2064.977363       0.0495
    LBFGS:   4  09:58:56    -2064.977446       0.0233
    LBFGS:   5  09:58:56    -2064.977460       0.0059
    $ ase gui H2O.traj@-1 -tg "a(1,0,2),d(0,1)"
    102.591620591 1.00793234388

.. highlight:: python

An example of creating an NWChem calculator in the python interface is::

  from ase.calculators.nwchem import NWChem

  calc = NWChem(label='calc/nwchem',
                dft=dict(maxiter=2000,
                         xc='B3LYP'),
                basis='6-31+G*')

Parameters
==========

.. highlight:: none

The NWChem calculator represents nested keyword blocks in the input file using
nested Python dictionaries. For example, consider the following block of input::

  memory 1024 mb

  dft
    xc B3LYP
    mult 2
    odft
    convergence energy 1e-5 density 1e-4 gradient 5e-3
  end

.. highlight:: python

This would be represented by the following keyword arguments to the NWChem
calculator::

  NWChem(memory='1024 mb',
         dft=dict(xc='B3LYP',
                  mult=2,
                  odft=None,
                  convergence=dict(energy=1e-5,
                                   density=1e-4,
                                   gradient=5e-3),
                  ),
         )


Most input files can be constructed in this way. The NWChem calculator also
has several special keywords which do not directly enter into the input file;
these are described in the table below

============== ======== =============== ==================================================
keyword        type     default value   description
============== ======== =============== ==================================================
``label``      ``str``  ``'nwchem'``    Used to name input and output files. Also,
                                        used as the default name of the perm and
                                        scratch directory, unless otherwise
                                        specified.
``theory``     ``str``  See description Theory specifies the kind of calculation
                                        you want to do. Currently supported
                                        values are ``dft``, ``scf``, ``mp2``,
                                        ``ccsd``, ``tce``, ``tddft``, ``pspw``,
                                        ``band``, and ``paw``. Other settings
                                        may work, but have not been tested.
                                        If not provided, the Calculator will
                                        attempt to guess based on the provided
                                        keywords.
``center``     ``bool`` ``False``       Whether NWChem should automatically
                                        center your atoms. You probably
                                        don't want to change this.
``autosym``    ``bool`` ``False``       Whether NWChem should automatically
                                        symmetrize your atoms. You probably
                                        don't want to change this.
``autoz``      ``bool`` ``False``       Whether NWChem should automatically
                                        generate a Z-matrix for your system.
``basis``      ``str``  ``3-21G``       If provided with a string, the
                                        specified basis set will be used for
                                        all atoms. Alternatively, you can set
                                        element-specific basis sets by passing
                                        a dict, e.g. ``basis=dict(C='6-31G', O='3-21G')``.
``basispar``   ``str``  ``''``          Additional keywords to go in the
                                        first line of the ``basis`` block.
``task``       ``str``  See description What kind of calculation is to be
                                        performed, e.g. ``'energy'``,
                                        ``'gradient'``, or ``'optimize'``.
                                        If not provided, it will be
                                        automatically determined by the
                                        Calculator.
``symmetry``   ``str``  ``''``          The symmetry group or number of your
                                        system.
``geompar``    ``str``  ``''``          Additional keywords to go in the first
                                        line of the ``geometry`` block.
``set``        ``dict`` ``dict()``      A set of keys and values to be added
                                        directly to the NWChem rtdb. This
                                        isn't necessary for most commonly
                                        done tasks, but it is required for
                                        certain functionality in plane-wave
                                        mode.
``restart_kw`` ``str``  ``'start'``     If set to ``'restart'``, then the calculation
                                        will attempt to restart from whatever
                                        restart files are available as
                                        specified in the .db file (may need to
                                        be combined with other keywords).
============== ======== =============== ==================================================

.. autoclass:: NWChem
