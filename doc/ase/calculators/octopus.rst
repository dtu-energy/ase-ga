.. module:: ase.calculators.octopus

=======
Octopus
=======


Introduction
============

Octopus_ is a density functional theory code focusing on
time-dependent simulations.  It supports several other calculation
modes as well.

This page documents the ASE interface to Octopus.  The interface is
written for Octopus tetricus (svn version 14450 or newer) but may be
compatible with older versions as well.

.. _Octopus: http://tddft.org/programs/octopus

Examples
========

Structure optimization
----------------------

Structure optimization of molecule:


.. literalinclude:: octopus_h2o.py

Numerical parameters can be specified as strings as well.

If the atoms have any cell vectors, then
``BoxShape`` will be ``parallelepiped`` unless user specifies
another ``BoxShape``.


Periodic system
---------------

Calculate density of states of silicon:

.. literalinclude:: octopus_silicon.py

Note how a *block* such as ``KPointsGrid`` is specified as a list of
lists.  In general, a block is a list of lists of strings.  Numerical
datatypes are converted to strings.

Time-dependent density functional theory
----------------------------------------

*TODO* Write script

Additional information
======================

The interface works by reading and writing files.  When triggering a
calculation, ASE writes an input file and proceeds to run Octopus.
Almost any keywords can be given to the Octopus calculator, and they
will be passed almost uncritically to Octopus.  This results in mostly
predictable but at times non-user-friendly behaviour.

ASE always works with ``Units = ev_angstrom``.  Keywords that attempt
to interfere with this will result in an error, or unspecified
behaviour if there are unhandled cases.

By default, Octopus output is written under the directory
``ink-pool``.  The ``label`` keyword can be used to specify a
different directory.

Hints
-----

Octopus input files can be loaded by the ASE GUI.  You can visualize
an Octopus input file by running ``ase gui inp``.

Bugs and misbehaviour
---------------------

Please report misbehaviour to the :ref:`ASE-developers <contact>`
unless you are certain that the behaviour is linked to Octopus and not
the ASE interface.

Below are listed some known bugs, issues, and limitations.  These may
or may not be fixed, depending on user response.

 * When parsing input files, arithmetic, calls to GNU GSL, or
   assignments from keywords are presently unsupported.  Only
   statements of the form ``keyword = value`` or blocks are supported.

 * Most Octopus keywords are passed directly to Octopus.  The ASE
   interface itself is not logically aware of their meaning.  Only
   those necessary to construct the Atoms are handled individually.
   There may exist keywords that affect the Atoms (the cell or
   geometry) which are not handled by the interface, and therefore
   result in confusing behaviour.  In particular, avoid changing the
   units.

 * Subsequent calculations always overwrite files.  This is probably
   fine for densities in a structure optimization (the density is
   reused on each step), but not useful for text output.  Therefore,
   avoid keywords like ``stdout='"out.log"'`` and redirect stdout and
   stderr by other means.

