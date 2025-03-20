.. module:: ase.filters

=======
Filters
=======

The Filter class
================

Constraints can also be applied via filters, which acts as a wrapper
around an atoms object. A typical use case will look like this::

   -------       --------       ----------
  |       |     |        |     |          |
  | Atoms |<----| Filter |<----| Dynamics |
  |       |     |        |     |          |
   -------       --------       ----------

and in Python this would be::

  >>> atoms = Atoms(...)
  >>> filter = Filter(atoms, ...)
  >>> dyn = Dynamics(filter, ...)


This class hides some of the atoms in an Atoms object.

.. class:: Filter(atoms, indices=None, mask=None)

You must supply either the indices of the atoms that should be kept
visible or a mask. The mask is a list of booleans, one for each atom,
being true if the atom should be kept visible.

Example of use::

  >>> from ase import Atoms
  >>> from ase.filters import Filter
  >>> atoms=Atoms(positions=[[ 0    , 0    , 0],
  ...                        [ 0.773, 0.600, 0],
  ...                        [-0.773, 0.600, 0]],
  ...             symbols='OH2')
  >>> f1 = Filter(atoms, indices=[1, 2])
  >>> f2 = Filter(atoms, mask=[0, 1, 1])
  >>> f3 = Filter(atoms, mask=[a.z == 1 for a in atoms])
  >>> f1.get_positions()
  [[ 0.773  0.6    0.   ]
   [-0.773  0.6    0.   ]]

In all three filters only the hydrogen atoms are made
visible.  When asking for the positions only the positions of the
hydrogen atoms are returned.


The UnitCellFilter class
========================

The unit cell filter is for optimizing positions and unit cell
simultaneously.  Note that :class:`ExpCellFilter` will probably
perform better.

.. autoclass:: UnitCellFilter

The StrainFilter class
======================

The strain filter is for optimizing the unit cell while keeping
scaled positions fixed.

.. autoclass:: StrainFilter


The ExpCellFilter class
=======================

The exponential cell filter is an improved :class:`UnitCellFilter`
which is parameter free.

.. autoclass:: ExpCellFilter

The FrechetCellFilter class
===========================

.. autoclass:: FrechetCellFilter
