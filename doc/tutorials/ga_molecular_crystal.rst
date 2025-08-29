.. _ga_molecular_crystal_tutorial:

=========================================================
Genetic algorithm search for molecular crystal structures
=========================================================

Similar to the :ref:`ga_bulk_tutorial` tutorial, this one illustrates how
to search for stable crystal structures using the GA, but this time for a
molecular crystal (solid dinitrogen, :mol:`N_2`). For such materials, the
global optimization problem can be strongly simplified by preserving
**molecular identity**. By limiting the search to crystals consisting of
:mol:`N_2` molecules, compared to individual :mol:`N` atoms (which may
or may not form :mol:`N_2` molecules), the number of degrees of freedom
is considerably reduced.

.. note::

  This assumption, of course, also needs to be checked --
  at high pressures, for example, nitrogen is
  :doi:`known to polymerize <10.1038/nmat1146>` and so
  will no longer consist of individual :mol:`N_2` molecules!

The same approach can also be used e.g. for compounds with polyatomic ions
or for adsorbates on surfaces. The main differences with the 'atomic' GA
therefore lie in the identity preservation during the various stages of the GA:

* The random initial structures need to be constructed from the appropriate
  building blocks (here: :mol:`N_2` molecules with a N-N bond length as in the
  gas phase).

* The genetic operators must leave the internal geometry of the polyatomic
  building blocks intact. In the present implementation, this is achieved
  by looking at the *tags* attribute of the parent Atoms objects:
  atoms belong to the same molecule if they have the same tag.
  These tags are either initialized by the initial structure generator, or
  can be added manually using the :meth:`Atoms.set_tags` method.

* The local optimization procedure may, if desired, change these internal
  geometries (i.e. the N-N bond lengths in our example). However, care may
  need to be taken that molecular identity is preserved also here (e.g. by
  applying certain constraints, see :mod:`ase.constraints`).


.. seealso::

  The :ref:`mhtutorial` tutorial addresses molecular identity
  in a different global optimization approach.


Initial population
==================

As usual, a GA run starts with generating a series of random initial
structures. The following script (:download:`ga_molecular_crystal_start.py`)
creates a :file:`gadb.db` database containing 10 structures, each
consisting of 8 :mol:`N_2` molecules.

.. literalinclude:: ga_molecular_crystal_start.py


Run the GA search
=================

Now the actual search can begin (:download:`ga_molecular_crystal_run.py`),
which should only take about five minutes to complete.

The ``relax`` function, which performs the variable-cell local optimization,
is imported from :download:`ga_molecular_crystal_relax.py`. For the purpose
of this tutorial, a simple interatomic potential is used where the
intramolecular interaction is described with a harmonic potential and the
intermolecular interactions with a Lennard-Jones potential.

.. note::

  Solid :mol:`N_2` takes on a `variety of crystal structures
  <https://en.wikipedia.org/wiki/Solid_nitrogen#Crystal_structure>`__.
  This short and simplified GA example is not intended as a
  thorough investigation of this variety. Typical structures produced
  by this example are e.g. cubic or hexagonal close packed with
  intermolecular distances on the order of 3.7 Å.

.. literalinclude:: ga_molecular_crystal_run.py
