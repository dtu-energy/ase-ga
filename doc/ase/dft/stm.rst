.. _stm:

STM images
==========

The STM is a revolutionary experimental surface probe that has
provided direct local insight into the surface electronic
structure. Sometimes the interpretation of STM topographs are not
straightforward and therefore theoretically modeled STM images may
resolve conflicting possibilities and point to an underlying atomistic
model. ASE includes python modules for generating
Tersoff-Hamann STM topographs.

The calculated tunneling current will be proportional to:

.. math::

    \int_{\epsilon_F}^{\epsilon_F+eV} \sum_{kn}
    w_{\mathbf k} |\Psi_{\mathbf k n}(\mathbf r)|^2
    \delta(\epsilon - \epsilon_{\mathbf k n}) d\epsilon,

where `V` is the bias voltage, `w_{\mathbf k}` is the `\mathbf k`-point weight
and `\Psi_{\mathbf k n}(\mathbf r)` is the wave function.

.. seealso::

   * `Tutorial using GPAW
     <https://gpaw.readthedocs.io/tutorialsexercises/electronic/stm/stm.html>`__

More details:

.. autoclass:: ase.dft.stm.STM
   :members:
