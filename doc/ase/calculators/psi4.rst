.. module:: ase.calculators.psi4

======
psi4
======

`psi4 <http://www.psicode.org/>`_ is an open source quatum chemistry code out of the Sherill Group at Georgia Tech. 

.. autoclass:: Psi4


Setup
=====

.. highlight:: bash

First we need to install psi4. There are `instructions <http://www.psicode.org/psi4manual/master/external.html>`_ available on their website for compiling the best possible version of psi4. However, the easiest way to obtain psi4 by obtaining the binary package from conda::

    conda install psi4 -c psi4; conda update psi4 -c psi4

The ase calculator operates using the psi4 python API, meaning that if psi4 is installed correctly you won't need to do anything else to get psi4 working. It is, however, recommended that you set up a psi4 scratch directory by setting the ``PSI_SCRATCH`` environment variable::

    export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files

This directory is where temporary electronic structure files will be written. It is important that this directory be located on the same machine as the calculation is being done to avoid slow read/write operations. This is set to ``/tmp`` by default. However, be aware that the ``/tmp`` directory might not be large enough.

Examples
========

.. highlight:: python

You can import psi4 and run it like any other calculator in ase::

    from ase.calculators.psi4 import Psi4
    from ase.build import molecule
    import numpy as np
    
    atoms = molecule('H2O')
    
    calc = Psi4(atoms = atoms,
            method = 'b3lyp',
            memory = '500MB' # this is the default, be aware!
            basis = '6-311g_d_p_')
    
    atoms.calc = calc
    print(atoms.get_potential_energy())
    print(atoms.get_forces())

However, once you have instantiated the psi4 ase calculator with an atoms object you can interact with the psi4 python API as well. The psi4 API is just an attribute of the psi4 ase calculator::

    calc.psi4.frequency('scf/cc-pvdz', molecule=calc.molecule, 
                        return_wfn=True, dertype=1)

This is not required though, as psi4 will act like any other ase calculator.

It should be noted that the ``method`` argument supports non-DFT methods (such as coupled cluster ``ccsd(t)``) as well. There is a great variety of `quatum <http://www.psicode.org/psi4manual/master/methods.html>`_ `methods <http://www.psicode.org/psi4manual/master/dft_byfunctional.html>`_ and `basis sets <http://www.psicode.org/psi4manual/master/basissets_tables.html>`_ to choose from.

Parallelization
===============

Psi4 runs on a single thread by default. However, you may increase the number of threads by passing in the ``num_threads`` argument, which can take either "max" or integer values.

Parameters
==========

The list of possible parameters and their defaults is shown below.
See the Psi4 documentation for full explanations of these different options.

================  ======== ======================== ============================
keyword           type     default value            description
================  ======== ======================== ============================
``label``         ``str``  ``'psi4-calc'``          Label for saved files.
``method``        ``str``  ``'hf'``                 Quantum Method or Functional
``charge``                 ``None``                 Charge
``basis``         ``str``  ``'aug-cc-pvtz'``        Basis set.
``memory``        ``str``  ``500MB``                The amount of memory allocated
                                                    to psi4
``num_threads``             ``1``                   The number of threads to run
                                                    psi4 on
``symmetry``      ``str``  ``'c1'``                 The symmetry of your system
``PSI_SCRATCH``   ``str``  ``/tmp``                 The scratch directory for
                                                    psi4
``multiplicity``  ``int``` ``None``                 The spin multiplicity of your
                                                    system

``reference``     ``str``  ``None``                 The reference wave function.
                                                    If you wish to run spin unrestricted
                                                    enter "uhf", otherwise, leave this
                                                    blank.
================  ======== ======================== ============================


