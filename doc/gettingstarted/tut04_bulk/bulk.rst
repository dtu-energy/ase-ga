Crystals and band structure
===========================

In this tutorial we calculate properties of crystals.

Setting up bulk structures
--------------------------

ASE provides three frameworks for setting up bulk structures:

 * :func:`ase.build.bulk`.  Knows lattice types
   and lattice constants for elemental bulk structures
   and a few compounds, but
   with limited customization.

 * :func:`ase.spacegroup.crystal`.  Creates atoms
   from typical crystallographic information such as spacegroup,
   lattice parameters, and basis.

 * :mod:`ase.lattice`.  Creates atoms explicitly from lattice and basis.

Let's run a simple bulk calculation.

.. admonition:: Exercise

   Use :func:`ase.build.bulk` to get a primitive cell of silver,
   then visualize it.

Silver is known to form an FCC structure, so presumably the function returned
a primitive FCC cell.  But it's always nice to be sure what we have in
front of us.  Can you recognize it as FCC?
You can e.g. use the ASE GUI to repeat the structure and recognize
the A-B-C stacking.
ASE should also be able to verify that it really is a primitive FCC cell
and tell us what lattice constant was chosen::

   print(atoms.cell.get_bravais_lattice())

Periodic structures in ASE are represented using ``atoms.cell`` and ``atoms.pbc``.
The cell is a :class:`~ase.cell.Cell` object which represents
the crystal lattice with three vectors.
``pbc`` is an array of three booleans indicating whether the
system is periodic in each direction.

Bulk calculation
----------------

For periodic DFT calculations we should generally use a number of
k-points which properly samples the Brillouin zone.
Many calculators including GPAW and Aims
accept the ``kpts`` keyword which can be a tuple such as
``(4, 4, 4)``.  In GPAW, the planewave mode
is very well suited for smaller periodic systems.
Using the planewave mode, we should also set a planewave cutoff (in eV)::

  from gpaw import GPAW, PW
  calc = GPAW(mode=PW(600), kpts=(8, 8, 8),
              setups={'Ag': '11'}, ...)

Here we have used the ``setups`` keyword to specify that we want the
11-electron PAW dataset instead of the default which has 17 electrons,
making the calculation faster.

(In principle, we should be sure to converge both kpoint sampling
and planewave cutoff -- I.e., write a loop and try different samplings
so we know both are good enough to accurately describe the quantity
we want.)


.. admonition:: Exercise

  Run a single-point calculation of bulk silver with GPAW.  Save the
  ground-state in GPAW's own format using ``calc.write('Ag.gpw')``.


Density of states
-----------------

Having saved the ground-state, we can reload it for ASE to extract
the density of states::

    import matplotlib.pyplot as plt
    from ase.dft.dos import DOS
    from gpaw import GPAW

    calc = GPAW('groundstate.gpw')
    dos = DOS(calc, npts=500, width=0)
    energies = dos.get_energies()
    weights = dos.get_dos()
    plt.plot(energies, weights)
    plt.show()

Calling the DOS class with ``width=0`` means ASE well calculate the DOS using
the linear tetrahedron interpolation method, which takes time but gives a nicer
representation.  We could also have given it
a nonzero width such as the default value of 0.1 (eV).  In that case it
would have used a simple Gaussian smearing with that width, but we would need
more k-points to get a plot of the same quality.

Note that the zero point of the energy axis is the Fermi energy.

.. admonition:: Exercise

   Plot the DOS.

   You probably recall that an Ag atom has 10 d electrons and one s electron.

   Which parts of the spectrum do you think
   originate (mostly) from s electrons?  And which parts (mostly) from
   d electrons?


Time for analysis.
As we probably know, the d-orbitals
in a transition metal atom are localized close to the nucleus while the
s-electron is much more delocalized.

In bulk systems, the s-states overlap a lot and therefore split into a
very broad band over a wide energy range.  d-states overlap much less
and therefore also split less:  They form a narrow band with a
very high DOS.  Very high indeed because there are 10 times as
many d electrons as there are s electrons.

So to answer the question, the d-band accounts for most of the states
forming the big, narrow chunk between -6.2 eV to -2.6 eV.  Anything outside
that interval is due to the much broader s band.

The DOS above the Fermi level may not be correct, since the SCF
convergence criterion (in this calculation)
only tracks the convergenece of occupied states.
Hence, the energies over the Fermi level 0 are probably wrong.


What characterises the noble metals Cu, Ag, and Au, is that the d-band
is fully occupied.  I.e.: The whole d-band lies below the Fermi level
(energy=0).
If we had calculated any other transition metal, the Fermi level would
lie somewhere within the d-band.

.. note::

   We could calculate the s, p, and d-projected DOS to see more
   conclusively which states have what character.
   In that case we should look up the GPAW documentation, or other
   calculator-specific documentation.  So let's not do that now.


Band structure
--------------

Let's calculate the band structure of silver.

First we need to set up a band path.  Our favourite image search
engine can show us some reference graphs.  We might find band
structures from both Exciting and GPAW with Brillouin-zone path
:math:`\mathrm{W L \Gamma X W K}`.  Luckily ASE knows these letters
and can also help us visualize the reciprocal cell::

  lat = atoms.cell.get_bravais_lattice()
  print(lat.description())
  lat.plot_bz(show=True)

In general, the :mod:`ase.lattice` module provides
:class:`~ase.lattice.BravaisLattice` classes used to represent each
of the 14 + 5 Bravais lattices in 3D and 2D, respectively.
These classes know
about the high-symmetry k-points and standard Brillouin-zone paths
(using the `AFlow <http://aflowlib.org/>`_ conventions).

.. admonition:: Exercise

   Build a band path for :math:`\mathrm{W L \Gamma X W K}`.
   You can use ``path = atoms.cell.bandpath(...)``
   --- see the :class:`~ase.cell.Cell` documentation for which
   parameters to supply.
   This gives us a :class:`~ase.dft.kpoints.BandPath` object.

You can ``print()`` the band path object to see some basic information about it,
or use its :meth:`~ase.dft.kpoints.BandPath.write` method to save the
band path to a json file such as :file:`path.json`.  Then visualize it
using the command::

  $ ase reciprocal path.json

Once we are sure we have a good path with a reasonable number of k-points,
we can run the band structure calculation.
How to trigger a band structure calculation depends
on which calculator we are using, so we would typically consult
the documentation for that calculator (ASE will one day provide
shortcuts to make this easier with common calculators)::

  calc = GPAW('groundstate.gpw')
  atoms = calc.get_atoms()
  path = atoms.cell.bandpath(<...>)
  calc = calc.fixed_density(kpts=path, symmetry='off')

We have here told GPAW to use our bandpath for k-points, not to
perform symmetry-reduction of the k-points, and to fix the electron
density.

Then we trigger a new calculation, which will be non-selfconsistent,
and extract and save the band structure::

  atoms.get_potential_energy()
  bs = calc.band_structure()
  bs.write('bs.json')

Again, the ASE command-line tool offers a helpful command to plot the
band structure from a file::

  $ ase band-structure bs.json

.. admonition:: Exercise

   Calculate, save, and plot the band structure of silver
   for the path :math:`\mathrm{W L \Gamma X W K}`.

You may need to zoom around a bit to see the whole thing at once.
The plot will show the Fermi level as a dotted line
(but does not define it as zero like the DOS plot before).
Looking at the band structure, we see the complex tangle of what must
be mostly d-states from before, as well as the few states with lower energy
(at the :math:`\Gamma` point) and higher energy (crossing the Fermi level)
attributed to s.


Equation of state
-----------------

We can find the optimal lattice parameter and calculate the bulk modulus
by doing an equation-of-state calculation.  This means sampling the energy
and lattice constant over a range of values to get the minimum as well
as the curvature, which gives us the bulk modulus.

The online ASE docs already provide a tutorial on how to do this, using
the empirical EMT potential:
https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html

.. admonition:: Exercise

   Run the EOS tutorial.


Complex crystals and cell optimisation
--------------------------------------

(If time is scarce, please consider skipping ahead to do the remaining
exercises before returning here.)

For the simple FCC structure we only have a single parameter, *a*, and
the EOS fit tells us everything there is to know.

For more complex structures we first need a more advanced framework
to build the atoms, such as the
:func:`ase.spacegroup.crystal` function.

The documentation helpfully tells us how to build a
rutile structure, saving us the trouble of looking up the atomic
basis and other crystallographic information.
Rutile is a common mineral form of :mol:`TiO_2`

.. admonition:: Exercise

   Build and visualize a rutile structure.

Let's uptimise the structure.
In addition to the positions, we must optimise the
unit cell which, being tetragonal, is characterised by the
two lengths *a* and *c*.

Optimising the cell requires the energy derivatives with respect to the
cell parameters accessible through the stress tensor.
``atoms.get_stress()`` calculates and returns the stress
as a vector of the 6 unique components (Voigt form).  Using it requires
that the attached calculator supports the stress tensor.  GPAW's planewave
mode does this.

The :class:`ase.constraints.ExpCellFilter` allows us to optimise
cell and positions simultaneously.  It does this by exposing the degrees
of freedom to the optimiser as if they were additional positions ---
hence acting as a kind of filter.
We use it by wrapping it around the atoms::

  from ase.optimize import BFGS
  from ase.constraints import ExpCellFilter
  opt = BFGS(ExpCellFilter(atoms), ...)
  opt.run(fmax=0.05)

.. admonition:: Exercise

  Use GPAW's planewave mode to optimize the rutile unit cell.
  You will probably need a planewave
  cutoff of at least 500 eV.  What are the optimised lattice constants
  *a* and *c*?

.. admonition:: Exercise

   Calculate the band structure of rutile.
   Does it agree with your favourite internet search engine?

Solutions
---------

Ag ground state:

.. literalinclude:: solution/bulk_part1_groundstate.py

Ag DOS:

.. literalinclude:: solution/bulk_part2_dos.py

Ag band structure:

.. literalinclude:: solution/bulk_part3_bandstructure.py

Rutile cell optimisation:

.. literalinclude:: solution/bulk_part4_cellopt.py


Rutile band structure:

.. literalinclude:: solution/bulk_part5_rutile_bands.py
