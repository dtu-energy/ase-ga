Nanoparticle
============

ASE provides a module, :mod:`ase.cluster`, to set up
metal nanoparticles with common crystal forms.
Please have a quick look at the documentation.

Build and optimise nanoparticle
-------------------------------

Consider :func:`ase.cluster.Octahedron`.  Aside from generating
strictly octahedral nanoparticles, it also offers a ``cutoff``
keyword to cut the corners of the
octahedron.  This produces "truncated octahedra", a well-known structural motif
in nanoparticles.  Also, the lattice will be consistent with the bulk
FCC structure of silver.

.. admonition:: Exercise

   Play around with :func:`ase.cluster.Octahedron` to produce truncated
   octahedra.  Set up a cuboctahedral
   silver nanoparticle with 55 atoms.  As always, verify with the ASE GUI that
   it is beautiful.

ASE provides a forcefield code based on effective medium theory,
:class:`ase.calculators.emt.EMT`, which works for the FCC metals (Cu, Ag, Au,
Pt, and friends).  This is much faster than DFT so let's use it to
optimise our cuboctahedron.

.. admonition:: Exercise

   Optimise the structure of our :mol:`Ag_55` cuboctahedron
   using the :class:`ase.calculators.emt.EMT`
   calculator.

Ground state
------------

One of the most interesting questions of metal nanoparticles is how
their electronic structure and other properties depend on size.
A small nanoparticle is like a molecule with just a few discrete energy
levels.  A large nanoparticle is like a bulk material with a continuous
density of states.  Let's calculate the Kohn--Sham spectrum (and density
of states) of our
nanoparticle.

As usual, we set a few parameters to save time since this is not a
real production calculation.
We want a smaller basis set
and also a PAW dataset with fewer electrons than normal.
We also want to use Fermi smearing since there could be multiple electronic
states near the Fermi level::

     from gpaw import GPAW, FermiDirac
     calc = GPAW(mode='lcao', basis='sz(dzp)', setups={'Ag': '11'},
                 occupations=FermiDirac(0,1))

These are GPAW-specific keywords --- with another code, those variables
would have other names.

.. admonition:: Exercise

   Run a single-point calculation of the optimised :mol:`Ag_55`
   structure with GPAW.

   After the calculation, dump the ground state to a file::

     calc.write('groundstate.gpw')


Density of states
-----------------

Once we have saved the ``.gpw`` file, we can write a new script
which loads it and gets the DOS::

  import matplotlib.pyplot as plt
  from gpaw import GPAW
  calc = GPAW('groundstate.gpw')
  energies, dos = calc.get_dos(npts=500, width=0.1)
  efermi = calc.get_fermi_level()

In this example, we sample the DOS using Gaussians of width 0.1 eV.
You will want to mark the Fermi level in the plot.  A good way
is to draw a vertical line: ``plt.axvline(efermi)``.

.. admonition:: Exercise

  Use matplotlib to plot the DOS as a function of energy, marking
  also the Fermi level.


.. admonition:: Exercise

   Looking at the plot, is this spectrum best understood as
   continuous or discrete?


The graph should show us that already with 55 atoms, the plentiful d
electrons are well on their way to forming a continuous band (recall
we are using 0.1 eV Gaussian smearing).  Meanwhile the energies of the
few s electrons split over a wider range, and we clearly see isolated
peaks: The s states are still clearly quantized and have significant
gaps.  What characterises the noble metals Cu, Ag, and Au,
is that their d band is fully occupied so that the Fermi level lies
among these s states.  Clusters with a
different number of electrons might have higher or lower Fermi level,
strongly affecting their reactivity.  We can conjecture that at 55
atoms, the properties of free-standing Ag nanoparticles are probably
strongly size dependent.

The above analysis is speculative.  To verify the analysis
we would want to calculate s, p, and d-projected DOS to see if our
assumptions were correct.  In case we want to go on doing this,
the GPAW documentation will be of help, see: `GPAW DOS <https://gpaw.readthedocs.io/tutorialsexercises/electronic/pdos/pdos.html#module-gpaw.dos>`__.


Solutions
---------

Optimise cuboctahedron:

.. literalinclude:: solution/Ag_part1_optimise.py

Calculate ground state:

.. literalinclude:: solution/Ag_part2_groundstate.py

Plot DOS:

.. literalinclude:: solution/Ag_part3_dos.py
