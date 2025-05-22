External calculators
--------------------

Many external calculators can be used with ASE, including GPAW_, Abinit_, Vasp_, Siesta_, `Quantum ESPRESSO`_, Asap_, LAMMPS_ and many more, see :ref:`supported calculators` for the full list.

.. _Asap: http://wiki.fysik.dtu.dk/asap
.. _GPAW: http://gpaw.readthedocs.io
.. _Siesta: http://www.icmab.es/siesta
.. _Abinit: https://www.abinit.org
.. _Vasp: https://www.vasp.at
.. _Quantum ESPRESSO: http://www.quantum-espresso.org/
.. _LAMMPS: http://lammps.sandia.gov/


Setting up an external calculator with ASE
==========================================

This tutorial will cover how to set up a basic calculation in ASE, using an external calculator.
We will be using the :class:`~ase.calculators.vasp.Vasp` class in this example, other calculators can be used in a similar manner. For details please refer to the documentation of the specific calculators :ref:`supported calculators`.

Important: ASE does not provide code or a license for VASP, and these must be acquired elsewhere.
ASE only creates an interface with VASP, so that you can use the ASE provided tools together with VASP.

Setting up
==========

The first step, is to tell ASE how to execute VASP, and where to find the pseudopotentials. You will need to have two environment variables defined:

.. highlight:: bash

::

   $ export ASE_VASP_COMMAND="mpiexec $HOME/vasp/bin/vasp_std"
   $ export VASP_PP_PATH=$HOME/vasp/mypps

The first environment variable :envvar:`ASE_VASP_COMMAND` is the default way to execute VASP, and should be defined in the same way, that you could normally execute a VASP run. Here we assume that the VASP executable, ``vasp_std``, is located in ``$HOME/vasp/bin``. Note, that if you want to execute VASP in parallel, this call should also include the MPI executable, which in this case is ``mpiexec``.

The second variable, :envvar:`VASP_PP_PATH`, is the path to the VASP pseudopotentials.

An additional (optional) variable for the :file:`vdw_kernel.bindat` file, which is required when doing van der Waals calculations, where ``luse_vdw=True``.

.. highlight:: bash

::

   $ export ASE_VASP_VDW=$HOME/<path-to-vdw_kernel.bindat-folder>

Note, that this should target the folder, and not the file itself.


Your first run
==============

Now that ASE knows how to execute VASP, we can try setting up a simple calculation. First we set up an atoms object

.. code-block:: python

    from ase.build import molecule

    atoms = molecule('N2')
    atoms.center(vacuum=5)

To perform a VASP DFT calculation, we now set up a calculator object.

.. code-block:: python

    from ase.calculators.vasp import Vasp

    calc = Vasp(xc='pbe',  # Select exchange-correlation functional
                encut=400, # Plane-wave cutoff
                kpts=(1, 1, 1)) # k-points

    atoms.calc = calc
    en = atoms.get_potential_energy()  # This call will start the calculation
    print('Potential energy: {:.2f} eV'.format(en))

Which results in the following output::

    Potential energy: -16.59 eV


The flow of how ASE interfaces with VASP, is that ASE handles writing the input files, which are required for the run, and then executes the :envvar:`ASE_VASP_COMMAND`, i.e. executes VASP.
Once the VASP run is complete, ASE then reads all of the relevant files, in this case the ``OUTCAR``, ``vasprun.xml`` and ``CONTCAR``, and stores properties in the calculator object.

For more information on the capabilities of the VASP calculators, please refer to :ref:`vasp-calculator`.
For other calculators, please refer to the :ref:`calculators` page.
