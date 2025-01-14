.. module:: ase.calculators.qmmm

QMMM
====

There are three QM/MM calculators native to ASE:

=========================  ===================
Explicit Interaction QMMM  :class:`EIQMMM`
Simple, subtrative QMMM    :class:`SimpleQMMM`
Force-based qmmm           :class:`ForceQMMM`
=========================  ===================


Explicit Interaction QMMM
-------------------------

In Explicit Interaction QMMM, the QM and MM regions
are explicitly coupled with an electrostatic interaction term.
This requires that the electrostatic potential from the classical charges of the
MM subsystem is fed into the QM calculator. This is built into GPAW_. More info
:doi:`in this paper <10.1021/acs.jctc.7b00621>`, which should be
cited if the method is used. 

Other ASE-calculators that currently support EIQMMM:

1. :mod:`DFTB+ <ase.calculators.dftb>`
2. :mod:`CRYSTAL <ase.calculators.crystal>`
3. :mod:`TURBOMOLE <ase.calculators.turbomole>`
4. :mod:`ORCA <ase.calculators.orca>` 

.. _GPAW: https://gpaw.readthedocs.io

.. seealso::

    The :ref:`qmmm` tutorial, on how to use the Explicit Interaction QMMM
    calculator

.. autoclass:: EIQMMM

Here, you need to specify the interaction::

    from ase.calculators.qmmm import EIQMMM, LJInteraction
    from ase.calculators.tip3p import epsilon0, sigma0
    lj = LJInteraction({'OO': (epsilon0, sigma0)})
    atoms.calc = EIQMMM([0, 1, 2],
                        QMCalculator(...),
                        MMCalculator(...),
                        interaction=lj)

For Lennard-Jones type of interactions you can use:

.. autoclass:: LJInteractions

Or, for couplings requiring more generality, you should use:

.. autoclass:: LJInteractionsGeneral

You can control how the QM part is embedded in the MM part by supplying your
own embedding object when you construct the :class:`EIQMMM` instance.  The
Embedding object will be specific to the QM calculator you want to use.  The
default is this one:

.. autoclass:: Embedding

Simple, subtractive QMMM calculations
-------------------------------------

This QM/MM calculator is similar to the original ONIOM model, doing
simple, subtractive QM/MM between any two calculators.

.. autoclass:: SimpleQMMM

This type of QMMM can combine any pair of ASE calculators::

    from ase.calculators.qmmm import SimpleQMMM
    atoms = ...
    atoms.calc = SimpleQMMM([0, 1, 2],
                            QMCalculator(...),
                            MMCalculator(...))

where ``[0, 1, 2]`` would select the first three atoms for the QM-part.


Force-based QM/MM
-----------------

This QM/MM calculator mixes forces from any pair of ASE calculators.
A finite buffer is added around the core QM region to ensure accurate forces; careful testing
of the required buffer size is required. See
:doi:`N. Bernstein, J. R. Kermode, and G. Csányi, Rep. Prog. Phys. 72, 026501 (2009) <10.1088/0034-4885/72/2/026501>`
for a review of force-based QM/MM approaches, which should be cited if this method is used,
and :doi:`T. D. Swinburne and J. R. Kermode, Phys. Rev. B 96, 144102 (2017) <10.1103/PhysRevB.96.144102>`
for an application which used this implementation.

.. autoclass:: ForceQMMM

Basic usage is as follows::

    from ase.calculators.qmmm import ForceQMMM
    atoms = ...
    qm_mask = ...
    atoms.calc = ForceQMMM(qm_mask,
                           QMCalculator(...),
                           MMCalculator(...),
                           buffer_width=...)

See :git:`ase/test/forcefields/test_forceqmmm.py` test-case for a complete example.
