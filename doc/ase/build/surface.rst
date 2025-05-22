========
Surfaces
========

.. currentmodule:: ase.build

.. _surfaces:

Common surfaces
===============

A number of utility functions are provided to set up
the most common surfaces, to add vacuum layers, and to add adsorbates
to a surface.  In general, all surfaces can be set up with
the modules described in the section :ref:`general-crystal-section`, but these
utility functions make common tasks easier.


Example
-------

To setup an Al(111) surface with a hydrogen atom adsorbed in an on-top
position::

  from ase.build import fcc111
  slab = fcc111('Al', size=(2,2,3), vacuum=10.0)

This will produce a slab 2x2x3 times the minimal possible size, with a
(111) surface in the z direction.  A 10 Å vacuum layer is added on
each side.

To set up the same surface with with a hydrogen atom adsorbed in an on-top
position 1.5 Å above the top layer::

  from ase.build import fcc111, add_adsorbate
  slab = fcc111('Al', size=(2,2,3))
  add_adsorbate(slab, 'H', 1.5, 'ontop')
  slab.center(vacuum=10.0, axis=2)

Note that in this case it is probably not meaningful to use the vacuum
keyword to fcc111, as we want to leave 10 Å of vacuum *after* the
adsorbate has been added. Instead, the :meth:`~ase.Atoms.center` method
of the :class:`~ase.Atoms` is used
to add the vacuum and center the system.

The atoms in the slab will have tags set to the layer number: First layer
atoms will have tag=1, second layer atoms will have tag=2, and so on.
Adsorbates get tag=0:

>>> print(atoms.get_tags())
[3 3 3 3 2 2 2 2 1 1 1 1 0]

This can be useful for setting up :mod:`ase.constraints` (see
:ref:`diffusion tutorial`).


Utility functions for setting up surfaces
-----------------------------------------

All the functions setting up surfaces take the same arguments.

*symbol*:
  The chemical symbol of the element to use.

*size*:
  A tuple giving the system size in units of the minimal unit cell.

*a*:
  (optional) The lattice constant.  If specified, it overrides the
  expermental lattice constant of the element.  Must be specified if
  setting up a crystal structure different from the one found in
  nature.

*c*:
  (optional) Extra HCP lattice constant.  If specified, it overrides the
  expermental lattice constant of the element.  Can be specified if
  setting up a crystal structure different from the one found in
  nature and an ideal `c/a` ratio is not wanted (`c/a=(8/3)^{1/2}`).

*vacuum*:
  The thickness of the vacuum layer.  The specified amount of
  vacuum appears on both sides of the slab.  Default value is None,
  meaning not to add any vacuum.  In that case the third axis perpendicular to
  the surface will be undefined (``[0, 0, 0]``) or left at its intrinsic
  bulk value if requested (see *periodic*).  Some calculators can work
  with undefined axes as long as the :attr:`~ase.Atoms.pbc` flag is set to
  ``False`` along that direction.

*orthogonal*:
  (optional, not supported by all functions). If specified and true,
  forces the creation of a unit cell with orthogonal basis vectors.
  If the default is such a unit cell, this argument is not supported.

*periodic*:
  (optional) Produce a bulk system.  Defaults to False.  If true, sets
  boundary conditions and cell constantly with the corresponding bulk
  structure.  Useful for stacking multiple different surfaces.  The
  system will be fully equivalent to the bulk material only if the
  number of layers is consistent with the crystal stacking.

Each function defines a number of standard adsorption sites that can
later be used when adding an adsorbate with
:func:`ase.build.add_adsorbate`.


The following functions are provided
````````````````````````````````````

.. autofunction:: fcc100
.. autofunction:: fcc110
.. autofunction:: bcc100
.. autofunction:: hcp10m10
.. autofunction:: diamond100

These always give orthorhombic cells:

==========  ============
fcc100      |fcc100|
fcc110      |fcc110|
bcc100      |bcc100|
hcp10m10    |hcp10m10|
diamond100  |diamond100|
==========  ============


.. autofunction:: fcc111
.. autofunction:: fcc211
.. autofunction:: bcc110
.. autofunction:: bcc111
.. autofunction:: hcp0001
.. autofunction:: diamond111

These can give both non-orthorhombic and orthorhombic cells:

===========  ===============  ===============
fcc111       |fcc111|         |fcc111o|
fcc211       not implemented  |fcc211o|
bcc110       |bcc110|         |bcc110o|
bcc111       |bcc111|         |bcc111o|
hcp0001      |hcp0001|        |hcp0001o|
diamond111   |diamond111|     not implemented
===========  ===============  ===============

The adsorption sites are marked with:

=======  ========  =====  =====  ========  ===========  ==========
ontop    hollow    fcc    hcp    bridge    shortbridge  longbridge
|ontop|  |hollow|  |fcc|  |hcp|  |bridge|  |bridge|     |bridge|
=======  ========  =====  =====  ========  ===========  ==========

.. |ontop|    image:: ontop-site.png
.. |hollow|   image:: hollow-site.png
.. |fcc|      image:: fcc-site.png
.. |hcp|      image:: hcp-site.png
.. |bridge|   image:: bridge-site.png
.. |fcc100|   image:: fcc100.png
.. |fcc110|   image:: fcc110.png
.. |bcc100|   image:: bcc100.png
.. |fcc111|   image:: fcc111.png
.. |bcc110|   image:: bcc110.png
.. |bcc111|   image:: bcc111.png
.. |hcp0001|  image:: hcp0001.png
.. |fcc111o|  image:: fcc111o.png
.. |fcc211o|  image:: fcc211o.png
.. |bcc110o|  image:: bcc110o.png
.. |bcc111o|  image:: bcc111o.png
.. |hcp0001o| image:: hcp0001o.png
.. |hcp10m10| image:: hcp10m10.png
.. |diamond100| image:: diamond100.png
.. |diamond111| image:: diamond111.png

.. autofunction:: mx2

.. image:: mx2.png

.. autofunction:: graphene

.. image:: graphene.png

Create root cuts of surfaces
````````````````````````````

To create some more complicated cuts of a standard surface, a root cell
generator has been created.  While it can be used for arbitrary cells,
some more common functions have been provided.

.. autofunction:: fcc111_root

.. autofunction:: hcp0001_root

.. autofunction:: bcc111_root

If you need to make a root cell for a different cell type, you can simply
supply a primitive cell of the correct height.  This primitive cell can be
any 2D surface whose normal points along the Z axis.  The cell's contents
can also vary, such as in the creation of an alloy or deformation.

.. autofunction:: ase.build.root_surface

The difficulty with using these functions is the requirement to know the
valid roots in advance, but a function has also been supplied to help with
this.  It is helpful to note that any primitive cell with the same cell
shape, such as the case with the fcc111 and bcc111 functions, will have the
same valid roots.

.. autofunction:: ase.build.root_surface_analysis

An example of using your own primitive cell::

  from ase.build import fcc111, root_surface
  atoms = fcc111('Ag', (1, 1, 3))
  atoms = root_surface(atoms, 27)

.. image:: fcc111_root.png


Adding adsorbates
-----------------

After a slab has been created, a vacuum layer can be added.  It is
also possible to add one or more adsorbates.

.. autofunction:: ase.build.add_adsorbate
.. autofunction:: ase.build.add_vacuum


.. _general-surface-section:

Create specific non-common surfaces
===================================

In addition to the most normal surfaces, a function has been
constructed to create more uncommon surfaces that one could be
interested in.  It is constructed upon the Miller Indices defining the
surface and can be used for both fcc, bcc and hcp structures.  The
theory behind the implementation can be found here:
:download:`general_surface.pdf`.

.. autofunction:: ase.build.surface


Example
-------

To setup a Au(211) surface with 9 layers and 10 Å of vacuum:

.. literalinclude:: general_surface.py
    :start-after: Au example
    :end-before: Mo example

This is the easy way, where you use the experimental lattice constant
for gold bulk structure.  You can write::

    from ase.visualize import view
    view(s1)

or simply ``s1.edit()`` if you want to see and rotate the structure.

.. image:: s1.png

Next example is a molybdenum bcc(321) surface where we decide what
lattice constant to use:

.. literalinclude:: general_surface.py
    :start-after: Mo example
    :end-before: Pt3Rh example

.. image:: s2.png

As the last example, creation of alloy surfaces is also very easily
carried out with this module.  In this example, two :mol:`Pt_3Rh`
fcc(211) surfaces will be created:

.. literalinclude:: general_surface.py
    :start-after: Pt3Rh example
    :end-before: Done

|s3| |s4|

.. |s3| image:: s3.png
.. |s4| image:: s4.png
