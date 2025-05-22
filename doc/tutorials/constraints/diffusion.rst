.. _constraints diffusion tutorial:

=======================================================
Surface diffusion energy barriers using ASE constraints
=======================================================

In this tutorial, we will calculate the energy barrier that was found
using the :mod:`NEB <ase.mep.neb>` method in the :ref:`diffusion tutorial`
tutorial.  Here, we use a simple :class:`~ase.constraints.FixedPlane`
constraint that forces the Au atom to relax in the *yz*-plane only:

.. literalinclude:: diffusion4.py

The result can be analysed with the command :command:`ase gui mep?.traj -n
-1` (choose :menuselection:`Tools --> NEB`).  The barrier is found to
be 0.35 eV - exactly as in the :ref:`NEB <diffusion tutorial>`
tutorial.

Here is a side-view of the path (unit cell repeated twice):

.. image:: diffusion-path.png


.. seealso::

   * :mod:`ase.mep.neb`
   * :mod:`ase.constraints`
   * :ref:`diffusion tutorial`
   * :func:`~ase.build.fcc100`
