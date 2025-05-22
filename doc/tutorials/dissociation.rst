.. _neb2:
.. _mep2:

===============================================
Dissociation of a molecule using the NEB method
===============================================

In this tutorial we provide an illustrative
example of a nudged-elastic band (NEB) calculation.
For more information on the NEB technique, see :mod:`ase.mep.neb`.
We consider the dissociation of a nitrogen molecule 
on the Cu (111) surface.

The first step is to find the relaxed structures
of the initial and final states.

.. literalinclude:: N2Cu-Dissociation1.py

Having obtained these structures we set up an NEB
calculation with 9 images.  Using :func:`~neb.interpolate()`
provides a guess for the path between the initial
and final states.  We perform the relaxation of the images
and obtain the intermediate steps.

.. literalinclude:: N2Cu-Dissociation2.py

After the calculation is complete, the energy difference
with respect to the initial state is given for each image,
as well as the distance between the N atoms.
