===================
Nudged elastic band
===================

.. module:: ase.mep.neb
   :synopsis: Nudged Elastic Band method.

The Nudged Elastic Band method is a technique for finding transition paths
(and corresponding energy barriers) between given initial and final states.
The method involves constructing a "chain" of "replicas" or "images" of the
system and relaxing them in a certain way.

Relevant literature References:


1. H. Jonsson, G. Mills, and K. W. Jacobsen, in 'Classical and Quantum
   Dynamics in Condensed Phase Systems', edited by B. J. Berne,
   G. Cicotti, and D. F. Coker, World Scientific, 1998 [standard
   formulation]

2. 'Improved Tangent Estimate in the NEB method for Finding Minimum
   Energy Paths and Saddle Points', G. Henkelman and H.
   Jonsson, J. Chem. Phys. 113, 9978 (2000) [improved tangent
   estimates]

3. 'A Climbing-Image NEB Method for Finding Saddle Points and Minimum
   Energy Paths', G. Henkelman, B. P. Uberuaga and H.
   Jonsson, J. Chem. Phys. 113, 9901 (2000)

4. 'Improved initial guess for minimum energy path calculations.',
   S. Smidstrup, A. Pedersen, K. Stokbro and H. Jonsson,
   J. Chem. Phys. 140, 214106 (2014)

5. 'Scaled and Dynamic Optimizations of Nudged Elastic Bands',
   P. Lindgren, G. Kastlunger and A. A. Peterson,
   J. Chem. Theory Comput. 15, 11, 5787-5793 (2019)


The NEB class
=============

This module defines one class:

.. autoclass:: NEB

Example of use, between initial and final state which have been previously
saved in A.traj and B.traj::

  from ase import io
  from ase.mep import NEB
  from ase.optimize import MDMin
  # Read initial and final states:
  initial = io.read('A.traj')
  final = io.read('B.traj')
  # Make a band consisting of 5 images:
  images = [initial]
  images += [initial.copy() for i in range(3)]
  images += [final]
  neb = NEB(images)
  # Interpolate linearly the potisions of the three middle images:
  neb.interpolate()
  # Set calculators:
  for image in images[1:4]:
      image.calc = MyCalculator(...)
  # Optimize:
  optimizer = MDMin(neb, trajectory='A2B.traj')
  optimizer.run(fmax=0.04)

Be sure to use the copy method (or similar) to create new instances
of atoms within the list of images fed to the NEB. Do *not* use something
like [initial for i in range(3)], as it will only create references to
the original atoms object.

Notice the use of the :meth:`~NEB.interpolate` method to obtain an
initial guess for the path from A to B.


Interpolation
=============

``NEB.interpolate()``

   Interpolate path linearly from initial to final state.

.. autofunction:: interpolate

   Interpolate path linearly from initial to final state. This standalone
   function can be used independently of the NEB class, but is functionally
   identical.

``NEB.interpolate(method='idpp')``

   Create an improved path from initial to final state using the IDPP approach
   [4]. This will start from an initial guess of a linear interpolation.

.. autofunction:: idpp_interpolate

   Generate an IDPP pathway from a set of images. This differs
   from above in that more IDPP-specific parameters can be specified,
   and an initial guess for the IDPP other than linear interpolation
   can be provided.

Only the internal images (not the endpoints) need have
calculators attached.


.. seealso::

   :mod:`ase.optimize`:
        Information about energy minimization (optimization). Note that you
        cannot use the default optimizer, BFGSLineSearch, with NEBs. (This is
        the optimizer imported when you import QuasiNewton.) If you would
        like a quasi-newton optimizer, use BFGS instead.

   :mod:`ase.calculators`:
        How to use calculators.

   :ref:`tutorials`:

        * :ref:`diffusion tutorial`
        * :ref:`neb2`
        * :ref:`idpp_tutorial`

.. note::

  If there are `M` images and each image has `N` atoms, then the NEB
  object behaves like one big Atoms object with `MN` atoms, so its
  :meth:`~ase.Atoms.get_positions` method will return a `MN \times 3`
  array.


Trajectories
============

The code::

  from ase.optimize import BFGS
  opt = BFGS(neb, trajectory='A2B.traj')

will write all images to one file.  The Trajectory object knows about
NEB calculations, so it will write `M` images with `N` atoms at every
iteration and not one big configuration containing `MN` atoms.

The result of the latest iteration can now be analysed with this
command: :command:`ase gui A2B.traj@-5:`.

For the example above, you can write the images to individual
trajectory files like this::

  for i in range(1, 4):
      opt.attach(io.Trajectory('A2B-%d.traj' % i, 'w', images[i]))

The result of the latest iteration can be analysed like this:

.. highlight:: bash

::

  $ ase gui A.traj A2B-?.traj B.traj -n -1

.. highlight:: python


Restarting
==========

Restart the calculation like this::

  images = io.read('A2B.traj@-5:')



Climbing image
==============

The "climbing image" variation involves designating a specific image to behave
differently to the rest of the chain: it feels no spring forces, and the
component of the potential force parallel to the chain is reversed, such that
it moves towards the saddle point. This depends on the adjacent images
providing a reasonably good approximation of the correct tangent at the
location of the climbing image; thus in general the climbing image is not
turned on until some iterations have been run without it (generally 20% to 50%
of the total number of iterations).

To use the climbing image NEB method, instantiate the NEB object like this::

  neb = NEB(images, climb=True)

.. note::

  Quasi-Newton methods, such as BFGS, are not well suited for climbing image
  NEB calculations. FIRE have been known to give good results, although
  convergence is slow.


Scaled and dynamic optimizations
================================

.. autoclass:: ase.mep.dyneb.DyNEB

The convergence of images is often non-uniform, and a large fraction of
computational resources can be spent calculating images that are below
the convergence criterion. This can be avoided with a dynamic optimization
method that monitors the convergence of each image. Dynamic optimization
is implemented as a subclass of the NEB method::

    from ase.mep import DyNEB
    neb = DyNEB(images, fmax=0.05, dynamic_relaxation=True)

where ``fmax`` must be identical to the ``fmax`` of the optimizer.

.. note::

    Dynamic optimization only works efficiently in series, and will not result
    in reduced computational time when resources are parallelized over images.
    ``dynamic_relaxation=False`` reverts to the default NEB implementation.

The saddle point is the important result of an NEB calculation, and the other
interior images are typically not used in subsequent analyses. The
convergence criteria can be scaled to focus on the saddle point and increase
the tolerance in other regions of the PES. Convergence scaling is implemented
as::

  neb = DyNEB(images, fmax=0.05, dynamic_relaxation=True, scale_fmax=1.)

where the convergence criterion of each image is scaled based on the position
of the image relative to the highest point in the band. The rate (slope) of
convergence scaling is controlled by the keyword ``scale_fmax``.

.. note::

  A low scaling factor (``scale_fmax=1-3``) is often enough to significantly
  reduce the number of force calls required for convergence.

Parallelization over images
===========================

Some calculators can parallelize over the images of a NEB calculation.
The script will have to be run with an MPI-enabled Python interpreter
like GPAW_'s gpaw-python_.  All images exist on all processors, but
only some of them have a calculator attached::

  from ase.parallel import world
  from ase.calculators.emt import EMT
  # Number of internal images:
  n = len(images) - 2
  j = world.rank * n // world.size
  for i, image in enumerate(images[1:-1]):
      if i == j:
          image.calc = EMT()

Create the NEB object with ``NEB(images, parallel=True)``.
For a complete example using GPAW_, see here_.

.. _GPAW: https://gpaw.readthedocs.io/
.. _gpaw-python: https://gpaw.readthedocs.io/documentation/parallel_runs/parallel_runs.html
.. _here: https://gpaw.readthedocs.io/tutorialsexercises/moleculardynamics/neb/neb.html


Using Shared Calculators
========================

Some calculators may parallelize well outside of ASE, or constructing them
consumes resources, for this purpose the user can allow the usage of shared
calculators with the ``allow_shared_calculator`` parameter of NEB.


.. _nebtools:

Analysis of output
==================

A class exists to help in automating the analysis of NEB jobs. See the
:ref:`Diffusion Tutorial <diffusion tutorial>` for some examples of its use.

.. autoclass:: NEBTools
   :members:

.. highlight:: bash

You can also make NEB plots that show the relaxation of your trajectory
directly from the command line; this will output the plots into a single PDF::

    $ ase nebplot neb.traj

You can find more help with::

    $ ase nebplot -h

.. highlight:: python

AutoNEB
=======

.. autoclass:: ase.mep.autoneb.AutoNEB
