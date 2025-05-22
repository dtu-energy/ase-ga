.. module:: ase.calculators.lammps

==================
LAMMPS Calculators
==================

LAMMPS_ (Large-scale Atomic/Molecular Massively Parallel Simulator) is
a classical molecular dynamics code.

    LAMMPS has potentials for soft materials (biomolecules, polymers) and solid-state materials (metals, semiconductors) and coarse-grained or mesoscopic systems.
    It can be used to model atoms or, more generically, as a parallel particle simulator at the atomic, meso, or continuum scale.

There are two calculators that interface to the LAMMPS molecular
dynamics code that can be used to solve an atoms model for energy,
atom forces and cell stresses. They are:

1. :class:`ase.calculators.lammpsrun.LAMMPS`
which interfaces to LAMMPS via writing a
controlling input file that is then run automatically through LAMMPS
and the results read back in. These results are currently limited to
total energy, atomic forces and cell stress.

2. :class:`ase.calculators.lammpslib.LAMMPSlib`
which uses the python interface that comes with LAMMPS,
loads the LAMMPS program as a python library. The LAMMPSlib calculator
then creates a '.lmp' object which is a running LAMMPS subroutine, so
further commands can be sent to this object and executed until it is
explicitly closed. Any additional variables calculated by LAMMPS can
also be extracted. Note however, any mistakes in the code sent to the
LAMMPS routine will cause python to terminate. Further information on the
python interface of LAMMPS can be found at lammpspy_link_. Note that it can be
very beneficial to compile lammps with C++ exceptions. Otherwise there will be
no error messages upon crashes.

It should not matter which code you use, but if you want access to
more of LAMMPS internal variables or to perform a more complicated
simulation then use LAMMPSlib. It is important to know which code you
are using because *when* you make an error in the LAMMPS code,
debugging the is difficult and different for both calculators.

Both of these interfaces are still experimental code and any
problems should be reported to the ASE developers mailing list.

.. _LAMMPS: https://lammps.sandia.gov/
.. _lammpspy_link: https://docs.lammps.org/Python_head.html

.. toctree::

    lammpsrun
    lammpslib
