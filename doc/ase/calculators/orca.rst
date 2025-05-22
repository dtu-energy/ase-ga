.. module:: ase.calculators.orca

======
ORCA
======

`ORCA <https://orcaforum.kofo.mpg.de/app.php/portal>`_ is a computational chemistry code
that can do SCF, (TD)DFT, semi-empirical potentials, MP2, CASSCF, Coupled Cluster
calculations, and more.


It is closed source, but free for academic users. Register on the forum to receive
a download link for the binaries, as well as access to the latest manual.


Many input examples are available at the
`ORCA Input Library <https://sites.google.com/site/orcainputlibrary>`_.


.. highlight:: none

The :class:`ORCA` ASE-interface is very simple. Two keywords are defined::

  orcasimpleinput: str
      What you'd put after the "!" in an orca input file.

  orcablocks: str
      What you'd put in the "% ... end"-blocks.


The ASE-calculator also works with the
:mod:`~ase.calculators.qmmm.EIQMMM`-calculator
for QM/MM simulations (see :mod:`~ase.calculators.qmmm` for
more info).

Setup and usage
===============

Orca can be configured using the configfile like other calculators.
If you need to override it for programmatic control of the orca
command, you can manually create an ``OrcaProfile``::

  from ase.calculators.orca import OrcaProfile

  profile = OrcaProfile(command='/full/path/to/my/orca')
  calc = ORCA(profile=profile)

ORCA decides which sub-processes to parallelize via MPI by itself, so you'll
almost always want a string in your ``orcablocks`` specifying the number of
cores for the simulation, e.g.::

  from ase.calculators.orca import ORCA

  calc = ORCA(profile=MyOrcaProfile,
              orcasimpleinput='B3LYP def2-TZVP'
              orcablocks='%pal nprocs 16 end')

for a B3LYP/def2-TZVP calculation on 16 cores.

Class Definition
================

.. autoclass:: ORCA
