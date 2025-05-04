.. _wannier tutorial:

============================================
Partly occupied Wannier Functions
============================================

This tutorial walks through building **partly occupied Wannier
functions** with the :mod:`ase.dft.wannier` module
and the `GPAW <https://wiki.fysik.dtu.dk/gpaw/>`_ electronic structure code. 
For more information on the details of the method and the implementation, see

   | K. S. Thygesen, L. B. Hansen, and K. W. Jacobsen
   | :doi:`Partly occupied Wannier functions: Construction and applications <https://doi.org/10.1103/PhysRevB.72.125119 >`
   | Phys. Rev. B 72, 125119, (2005)

.. contents:: **Outline**
   :depth: 2
   :local:


Benzene molecule
================

Step 1 – Ground-state calculation
---------------------------------

Run the script below to obtain the ground-state density and the
Kohn–Sham (KS) orbitals. The result is stored in :file:`benzene.gpw`.

.. literalinclude:: benzene.py
   :language: python

Step 2 – Maximally localized WFs for the occupied subspace (15 WFs)
-------------------------------------------------------------------

There are 15 occupied bands in the benzene molecule. We construct one Wannier function per occupied band by setting
``nwannier = 15``. 
By calling ``wan.localize()``, the code attempts to minimize the spread functional using a gradient-descent algorithm. 
The resulting WFs are written to .cube files, which allows them to be inspected using e.g. VESTA.

.. literalinclude:: wannier_benzene.py
   :language: python

Step 3 – Adding three extra degrees of freedom (18 WFs)
-------------------------------------------------------

To improve localization we augment the basis with three extra Wannier functions - so-called *extra degrees of freedom*
(``nwannier = 18``, ``fixedstates = 15``). This will allow the Wannierization procedure to use the unoccupied states to minimize spread functional.

.. literalinclude:: wannier_benzene_with_edf.py
   :language: python


Step 4 – Spectral-weight analysis
---------------------------------

The script below projects the WFs on the KS eigenstates. You should see
the 15 lowest bands perfectly reconstructed (weight ≃ 1.0) while higher
bands are only partially represented.

.. literalinclude:: plot_spectral_weight.py
   :language: python


Polyacetylene chain (1-D periodic)
==================================

We now want to construct partially occupied Wannier functions to describe a polyacetylene chain.

Step 1 – Structure & ground-state calculation
---------------------------------------------

Polyacetylene is modelled as an infinite chain; we therefore enable
periodic boundary conditions along *x*.

.. literalinclude:: polyacetylene.py
   :language: python

Step 2 – Wannierization
-----------------------

We repeat the localization procedure, keeping the five lowest
bands fixed and adding one extra degree of freedom to aid localization.

.. literalinclude:: wannier_polyacetylene.py
   :language: python

Step 3 – High-resolution band structure
---------------------------------------

Using the Wannier Hamiltonian we can interpolate the band structure on a
fine 100-point *k* mesh and compare it to the original DFT result. 

.. literalinclude:: plot_band_structure.py
   :language: python

Within the fixed-energy window—that is, for energies below the fixed-energy line—the Wannier-interpolated bands coincide perfectly with the DFT reference (red circles). Above this window the match is lost, because the degrees of freedom deliberately mix several Kohn–Sham states to achieve maximal real-space localisation.