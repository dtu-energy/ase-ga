Molecular vibrations
====================

Let's calculate the vibrational modes of :mol:`H_2O`.

Consider the molecule at its equilibrium positions.  If we displace the
atoms slightly, the energy :math:`E` will increase, and restoring
forces will make the atoms oscillate in some pattern around the equilibrium.

We can Taylor expand the energy with respect to the 9 coordinates (generally :math:`3N` coordinates for a molecule with :math:`N` atoms), :math:`u_i`:

.. math::

   E = E_0 + \frac{1}{2}\sum_{i}^{3N} \sum_{j}^{3N} \frac{\partial^2 E}{\partial u_{i}\partial u_{j}}\bigg\rvert_0 (u_i - u_{i0}) (u_j - u_{j0}) + \cdots

Since we are expanding around the equilibrium positions, the energy
should be stationary and we can omit linear contributions.

The matrix of all the second derivatives is called the Hessian, :math:`\mathbf H`, and it expresses a linear system of differential equations

.. math::

  \mathbf{Hu}_k = \omega_k^2\mathbf{Mu}_k

for the vibrational eigenmodes :math:`u_k` and their frequencies
:math:`\omega_k` that will characterise the collective movement of the atoms.
In short, we need the eigenvalues and eigenvectors of the Hessian.

The elements of the Hessian can be approximated as

.. math::

  H_{ij} = \frac{\partial^2 E}{\partial u_{i}\partial u_{j}}\bigg\rvert_0 = -\frac{\partial F_{j}}{\partial u_{i}},

where :math:`F_j` are the forces.  Hence we calculate the derivative
of the forces using finite differences.
We need to displace each atom back and forth along each Cartesian direction,
calculating forces at each configuration to establish
:math:`H_{ij} \approx \Delta F_{j} / \Delta u_{i}`,
then get eigenvalues and vectors of that.


ASE provides the :class:`~ase.vibrations.Vibrations` class for this
purpose.  Note how the linked documentation contains an example for
the :mol:`N_2` molecule, which means we almost don't have to do any
work ourselves.  We just scavenge the online ASE
documentation like we always do, then hack as necessary until the thing runs.


.. admonition:: Exercise

   Calculate the vibrational frequencies of :mol:`H_2O` using GPAW in
   LCAO mode, saving the modes to trajectory files.  What are the
   frequencies, and what do the eigenmodes look like?

Since there are nine coordinates, we get nine eigenvalues and
corresponding modes.  However the three translational and three
rotational degrees of freedom will contribute six "modes" that do not
correspond to true vibrations.  In principle there are no restoring
forces if we translate or rotate the molecule, but these will
nevertheless have different energies (often imaginary) because of
various artifacts of the simulation such as the grid used to represent
the density, or effects of the simulation box size.

A solution and other comments to this exercise can be found on the
GPAW web page:

https://gpaw.readthedocs.io/tutorialsexercises/vibrational/vibrations/vibrations.html