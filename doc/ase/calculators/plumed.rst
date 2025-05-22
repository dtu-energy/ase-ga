.. module:: ase.calculators.plumed

======
PLUMED
======

.. image:: ../../static/plumed.png

Introduction
============

Plumed_ is an open source library that allows implementing several
kinds of enhanced sampling methods and contains a variety of tools to
analyze data obtained from molecular dynamics simulations. With this 
calculator, You can carry out biased simulations including metadynamics,
well-tempered metadynamics, among others. Besides, it is possible to 
compute a large set of collective variables that plumed 
has already implemented for being calculated on-the-fly in MD simulations 
or for postprocessing tasks.

.. _Plumed: https://www.plumed.org/ 

Installation
============
The ASE-Plumed calculator uses the python wrap of Plumed. An easy way to
install it using conda::

    conda install -c conda-forge py-plumed

However, the installation preferences could be easier to modify using
any of the others options presented in
`this page <https://www.plumed.org/doc-v2.7/user-doc/html/_installation.html#installingpython>`_.

Test the installation of plumed doing this:

    >>> from plumed import Plumed
    >>> Plumed()


Set-up
======

Typically, Plumed simulations need an external file, commonly called plumed.dat
for setting up its functions. In this ASE calculator interface, Plumed
information is given to the calculator through a list of strings containing the
lines that would be included in the plumed.dat file. For example::

    setup = [f"UNITS LENGTH=A TIME={1/(1000 * units.fs)} ENERGY={units.mol/units.kJ}",
             "d: DISTANCE ATOMS=1,2",
             "PRINT ARG=d STRIDE=10 FILE=COLVAR"]

For a complete explanation of the plumed keywords, visit
`the official web of plumed plug-in <https://www.plumed.org/doc>`_.


Units
"""""

Note that the first line of setup list of the previous example is referred to
units. That is because Plumed will consider all quantities of input and outputs in
`plumed internal units <https://www.plumed.org/doc-v2.8/user-doc/html/_u_n_i_t_s.html>`_.
Then, it is necessary to add this line in order to remain the units same as ASE.
You can ignore this line, but be aware of the units changes.


.. seealso::

    Visit 
    `this Metadynamics tutorial <https://gitlab.com/Sucerquia/ase-plumed_tutorial>`_
    for further explanation of the Plumed calculator.

Plumed Calculator Class
=======================

.. autoclass:: ase.calculators.plumed.Plumed

.. note::
   Periodic Boundary Conditions (PBC) fixed in ASE has nothing to do with PBC
   of PLUMED. If you set a cell in ASE, PLUMED will assume PBC in all
   directions -at least you specify something different in your plumed set up-
   independently of your :mod:`~ase.Atoms.set_pbc` election.
