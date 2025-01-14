.. _tutorials:

Tutorials
=========

Python
------

If you are not familiar with Python please read :ref:`what is python`.

.. toctree::
   :hidden:

   ../python

If your ASE scripts make extensive use of matrices you may want to familiarize yourself with :ref:`numpy`.

.. toctree::
   :hidden:

   ../numpy

ASE
---

Most of the tutorials will use the :mod:`EMT <ase.calculators.emt>` potential,
but any other :mod:`Calculator <ase.calculators>` could be plugged in instead.

Basic property calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   atomization
   eos/eos
   lattice_constant

Surface adsorption
^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   db/db

Global optimization
^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   minimahopping/minimahopping
   ga/ga_optimize
   ga/ga_fcc_alloys
   ga/ga_convex_hull
   ga/ga_bulk
   ga/ga_molecular_crystal

Calculating diffusion/dissociation properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   neb/diffusion
   constraints/diffusion
   dissociation
   neb/idpp
   selfdiffusion/al110

ASE database
^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   tut06_database/database

Molecular Dynamics
^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   md/md
   tipnp_equil/tipnp_equil
   acn_equil/acn_equil

Uncategorized
^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   defects/defects
   qmmm/qmmm
   dimensionality/dimensionality
   deltacodesdft/deltacodesdft
   wannier/wannier
   tut03_vibrations/vibrations


Further reading
---------------

For more details:

* Look at the documentation for the individual :ref:`modules <ase>`.
* Browse the :git:`source code <>` online.


Videos
------

The following video tutorials are available:

 - **Overview and installation of ASE**, by Anthony Goodrow (duration: ~5min 30sec; size: 26 MB) - en: |oi_en|

.. |oi_en| image:: ../static/United_States_of_America.png
   :target: https://wiki.fysik.dtu.dk/ase-files/oi_en.avi

.. |oi_cn| image:: ../static/China.png
   :target: https://wiki.fysik.dtu.dk/ase-files/oi_ch.avi

