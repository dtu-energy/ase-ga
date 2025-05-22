.. module:: ase.calculators.dftd3

=======
DFT-D3
=======


Introduction
============

The DFTD3_ calculator class wraps the 'dftd3' command line utility by
the research group of Stefan Grimme. This can be used to calculate classical
vdW dispersion corrections to a large number of common DFT functionals. This
calculator can be used in conjunction with other DFT calculators such as
GPAW to allow seamless calculation of dispersion-corrected DFT energies,
forces, and stresses.

.. _DFTD3: https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/

This is a list of all supported keywords and settings:

===========  =============  ==================================================
Keyword      Default value  Description
===========  =============  ==================================================
``xc``       ``'pbe'``      Use parameters optimized for the selected XC 
                            functional.
``func``     ``None``       Alternative to ``xc``. Use one or the other.
``grad``     ``True``       Enable or disable calculation of gradients
                            (forces, stress tensor).
``abc``      ``False``      Enable three-body ATM correction.
``cnthr``    40 Bohr        Cutoff radius for coordination number and
                            three-body calculations.
``cutoff``   95 Bohr        Cutoff radius for two-body dispersion
                            calculations.
``old``      ``False``      Enable older DFT-D2 dispersion correction method.
``damping``  ``'zero'``     Damping method. Valid options are ``'zero'``,
                            ``'bj'``, ``'zerom'``, and ``'bjm'``.
``tz``       ``False``      Custom parameters optimized for
                            triple-zeta basis sets.
``s6``                      Custom damping parameter used in all damping
                            methods.
``sr6``                     Custom damping parameter used in ``'zero'`` and
                            ``'zerom'`` damping methods.
``s8``                      Custom damping parameter used in all damping
                            methods.
``sr8``                     Custom damping parameter used in ``'zero'`` and
                            ``'zerom'`` damping methods.
``alpha6``                  Custom damping parameter used in all damping
                            methods.
``a1``                      Custom damping parameter used in ``'bj'`` and
                            ``'bjm'`` damping methods.
``a2``                      Custom damping parameter used in ``'bj'`` method.
``beta``                    Custom damping parameter used in ``'bjm'``
                            method.
===========  =============  ==================================================

Examples
========

Before using DFT-D3 be sure you have installed its current version (link_).

.. _link: https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/

DFTD3 can be used by itself to calculate only the vdW correction to a
system's energy, forces, and stress. Note that you should not use these
properties alone to perform dynamics, as DFTD3 is not a full classical
potential.

.. literalinclude:: dftd3_alone.py

If used in conjunction with a DFT calculator, DFTD3 returns
dispersion-corrected energies, forces, and stresses which can be used to
perform dynamics.

.. literalinclude:: dftd3_gpaw.py

Additional information
======================

This calculator works by writing either an ``xyz`` file (for non-periodic
systems) or a ``POSCAR`` file (for periodic systems), calling the
``dftd3`` executable, and parsing the output files created. It has been
written such that its interface should match that of the ``dftd3`` utility
itself as closely as possible, while minimizing the possibility of setting
redundant and contradictory options. For example, you can only select one
damping method, and the interface will sanity-check any provided custom
damping parameters.

Without any arguments, the DFTD3 will default to calculating the PBE-D3
dispersion correction with ``'zero'`` damping. If a DFT calculator is
attached, DFTD3 will attempt to glean the XC functional from the DFT
calculator. This will occasionally fail, as ``dftd3`` is very particular
about how the names of XC functionals are to be formatted, so in general
you should supply the XC functional to both the DFT calculator and the DFTD3
calculator.

Caveats
-------

The ``dftd3`` does not handle systems with only 1D- or 2D-periodic boundary
conditions. If your system has 1D or 2D PBC, DFTD3 will calculate the
dispersion correction as though it was fully 3D periodic.

If your system is very large, the dispersion calculation can take quite long,
especially if you are including three-body corrections (``abc=True``). For
highly parallel calculations, this may result in the dispersion correction
taking longer than the DFT calculation! This is because the ``dftd3`` utility
is not parallelized and will always run on a single core. Be sure to
benchmark this calculator interface on your system before deploying large,
heavily parallel calculations with it!
