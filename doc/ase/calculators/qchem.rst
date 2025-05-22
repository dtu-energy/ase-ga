.. module:: ase.calculators.qchem

======
Q-Chem
======

.. image:: ../../static/qchem.png
   :target: https://www.q-chem.com/
   :scale: 40

`Q-Chem <https://www.q-chem.com/>`_ is a comprehensive ab initio quantum
chemistry package for accurate predictions of molecular structures,
reactivities, and vibrational, electronic and NMR spectra.


Setup
=====

.. highlight:: bash

You first need to install a working copy of Q-Chem for ASE to call;
follow the instructions on the `Q-Chem website <https://www.q-chem.com/>`_.

The default command that ASE will use to start Q-Chem is
``qchem PREFIX.inp PREFIX.out``.


Examples
========

.. highlight:: python

An simple example of running a geometry optimization using the QChem calculator
in the python interface::

  from ase.build import molecule
  from ase.calculators.qchem import QChem
  from ase.optimize import LBFGS

  mol = molecule('C2H6')
  calc = QChem(label='calc/ethane',
               method='B3LYP',
               basis='6-31+G*')
  mol.calc = calc
  opt = LBFGS(mol)
  opt.run()

See the source code link below for further details.

.. autoclass:: QChem
