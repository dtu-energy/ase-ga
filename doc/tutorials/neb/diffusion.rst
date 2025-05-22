.. _diffusion tutorial:

=============================================================================
Surface diffusion energy barriers using the Nudged Elastic Band  (NEB) method
=============================================================================

First, set up the initial and final states:

|initial| |final|

.. literalinclude:: diffusion1.py

.. note::  Notice how the tags are used to select the constrained atoms

Now, do the NEB calculation:

.. literalinclude:: diffusion2.py

Visualize the results with::

    $ ase gui neb.traj@-5:

and select Tools->NEB.

|ts| |barrier|

You can also create a series of plots like above, that show the progression
of the NEB relaxation, directly at the command line::

    $ ase nebplot --share-x --share-y neb.traj

For more customizable analysis of the output of many NEB jobs, you can use
the :class:`ase.mep.NEBTools` class. Some examples of its use are below; the
final example was used to make the figure you see above.

.. literalinclude:: diffusion5.py

.. note::

   For this reaction, the reaction coordinate is very simple: The
   *x*-coordinate of the Au atom.  In such cases, the NEB method is
   overkill, and a simple constraint method should be used like in this
   tutorial: :ref:`constraints diffusion tutorial`.

.. seealso::

   * :mod:`ase.mep.neb`
   * :mod:`ase.constraints`
   * :ref:`constraints diffusion tutorial`
   * :func:`~ase.build.fcc100`



.. |initial| image:: diffusion-I.png
.. |final| image:: diffusion-F.png
.. |ts| image:: diffusion-T.png
.. |barrier| image:: diffusion-barrier.png


Restarting NEB
==============

Restart NEB from the trajectory file:

.. literalinclude:: diffusion4.py


Parallelizing over images with MPI
==================================

Instead of having one process do the calculations for all three
internal images in turn, it will be faster to have three processes do
one image each. In order to be able to run python with MPI
you need a special parallel python interpreter, for example gpaw-python.

The example below can then be run
with ``mpiexec -np 3 gpaw-python diffusion3.py``:

.. literalinclude:: diffusion3.py
