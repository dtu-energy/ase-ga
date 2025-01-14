.. module::  ase.calculators.emt
   :synopsis: Effective Medium Theory

==========================
Pure Python EMT calculator
==========================

The EMT potential is included in the ASE package in order to have a
simple calculator that can be used for quick demonstrations and
tests.

.. warning::

   If you want to do a real application using EMT, you should use the
   *much* more efficient implementation in the ASAP_ calculator.

.. autoclass:: ase.calculators.emt.EMT
   :class-doc-from: class

.. _ASAP: http://wiki.fysik.dtu.dk/asap

Theory
------

In the following, the seven parameters :math:`E_{0i}`, :math:`s_{0i}`,
:math:`V_{0i}`, :math:`\eta_{2i}`, :math:`\kappa_{i}`, :math:`\lambda_{i}`,
and :math:`n_{0i}` are specific for the species of atom :math:`i`.

Energy
^^^^^^

In the effective-medium theory (EMT), the energy is given by

.. math::
   E = \sum_{i=1}^{N} (E_{\mathrm{c},i} + E_{\mathrm{AS},i})

The cohesive function :math:`E_{\mathrm{c},i}` describes the energy in the
reference system, where we assume the face-centered cubic (fcc) structure and
given by

.. math::
   E_{\mathrm{c},i}
   = E_{0i} f(\lambda_i (s_{i} - s_{0i}))
   = E_{0i} f(\lambda_i \dot{s}_{i})

.. math::
   f(x) = (1 + x) \exp(-x)

where :math:`E_{0i}` is the cohesive energy, :math:`s_{0i}` is the Wigner–Seitz
radius in the equilibrium fcc state, and :math:`\lambda_i` is related to the
curvature of the energy–volume curve and thus to the bulk modulus.
:math:`s_i` is the neutral-sphere radius, and

.. math::
   \dot{s}_i
   = s_{i} - s_{0i}
   = - \frac{1}{\beta \eta_{2i}} \log \frac{\sigma_{1i}}{12 \gamma_{1i}}

where :math:`\beta` is the constant related to the Wigner–Seitz radius and the
first nearest neighbor distance (cf. `Tips`_).
:math:`\sigma_{1i}`` is given by

.. math::
   \sigma_{1i}
   = \sum_{j} \chi_{ij} w(r_{ij})
   \exp(- \eta_{2j} (r_{ij} - \beta s_{0j}))
   = \sum_{j} \dot{\sigma}_{1ij}^\mathrm{s}

The summation is over the neighbors of atom :math:`i`.
:math:`r_{ij}` is the distance of atoms :math:`i` and :math:`j` and given
using their position vectors :math:`\mathbf{r}_i` and :math:`\mathbf{r}_j` by

.. math::
   r_{ij} = |\mathbf{r}_{ij}| = |\mathbf{r}_{j} - \mathbf{r}_i|

:math:`\chi_{ij}` is given by

.. math::
   \chi_{ij} = \frac{n_{0j}}{n_{0i}}

The contribution from atom :math:`j` is given by

.. math::
   \dot{\sigma}_{1ij}^\mathrm{s} = \chi_{ij} w(r_{ij})
   \exp(- \eta_{2j} (r_{ij} - \beta s_{0j}))

For later convenience in `Forces`_, the contribution from atom :math:`i` to
atom :math:`j` is also written as;

.. math::
   \dot{\sigma}_{1ij}^\mathrm{o} = \chi_{ji} w(r_{ij})
   \exp(- \eta_{2i} (r_{ij} - \beta s_{0i}))

:math:`w(r)` is the smooth cutoff function given by

.. math::
   w(r) = \frac{1}{1 + \exp(a (r - r_\mathrm{c}))}

:math:`\gamma_{1i}` is a correction factor when considering beyond the first
nearest neighbor sites and given by (cf. `Tips`_)

.. math::
   \gamma_{1i} = \frac{1}{12} (
      n^\mathrm{1NN} w(d_0^\mathrm{1NN}) \exp(\eta_{2i} (d_0^\mathrm{1NN} - \beta s_{0i})) +
      n^\mathrm{2NN} w(d_0^\mathrm{2NN}) \exp(\eta_{2i} (d_0^\mathrm{2NN} - \beta s_{0i})) +
      n^\mathrm{3NN} w(d_0^\mathrm{3NN}) \exp(\eta_{2i} (d_0^\mathrm{3NN} - \beta s_{0i})) +
      \cdots
   )

which is :math:`1` when considering only up to the first nearest neighbors of
the equilibrium fcc structure.

The atomic-sphere correction :math:`E_{\mathrm{AS},i}` describes the derivation
from the reference fcc system and given by

.. math::
   E_{\mathrm{AS},i} = E_{\mathrm{AS},i}^{1} + E_{\mathrm{AS},i}^{2}

:math:`E_{\mathrm{AS},i}^{1}` is the pair interactions of the real system,
and :math:`E_{\mathrm{AS},i}^{2}` is the negative of the pair interactions of
the reference unary perfect fcc structure.
Both terms are described using the following pair interaction function;

.. math::
   V_{ij} (r) = -V_{0i} \cdot \frac{1}{\gamma_{2i}} \chi_{ij} w(r)
   \exp(-\frac{\kappa_{j}}{\beta}(r - \beta s_{0j}))

:math:`\gamma_{2i}` is a correction factor when considering beyond the first
nearest neighbor sites and given by (cf. `Tips`_)

.. math::
   \gamma_{2i} = \frac{1}{12} (
      n_\mathrm{1NN} w(d_\mathrm{1NN}) \exp(-\frac{\kappa_{2i}}{\beta} (r_\mathrm{1NN} - \beta s_{0i})) +
      n_\mathrm{2NN} w(d_\mathrm{2NN}) \exp(-\frac{\kappa_{2i}}{\beta} (r_\mathrm{2NN} - \beta s_{0i})) +
      n_\mathrm{3NN} w(d_\mathrm{3NN}) \exp(-\frac{\kappa_{2i}}{\beta} (r_\mathrm{3NN} - \beta s_{0i})) +
      \cdots
   )

Here, if we consider only up to the first nearest neighbors,

.. math::
   \gamma_{2i} \rightarrow 1

For :math:`E_{\mathrm{AS},i}^{2}`, only the interactions up to the first
nearest neighbors are considered, i.e., :math:`j = i` and
:math:`r_{ij} = d^\mathrm{1NN} = \beta s_{i}`. Thus,

.. math::
   E_{\mathrm{AS},i}^{2} = \frac{1}{2} n^\mathrm{1NN} V_{ii}(d^\mathrm{1NN})
   = -\frac{12}{2} V_{ii} (\beta s_{i})
   = 6 V_{0i} \exp(-\kappa_{i} \dot{s}_i)

The first term :math:`E_{\mathrm{AS},i}^{1}` is the pair interactions of the
real system. Here we consider the interactions up to a certain cutoff radius,
and we average the contribution from atom :math:`i` to atom :math:`j` and that
from atom :math:`j` to atom :math:`i`. Thus,

.. math::
   E_{\mathrm{AS},i}^{1} = \frac{1}{2}
   \sum_{j} \frac{1}{2} \left(V_{ij}(r_{ij}) + V_{ji}(r_{ij})\right)
   = - \frac{V_{0i}}{2 \gamma_{2i}} \cdot \frac{1}{2} \sum_{j}
   (\dot{\sigma}_{2ij}^\mathrm{s} + \dot{\sigma}_{2ij}^\mathrm{o})

where

.. math::
   \dot{\sigma}_{2ij}^\mathrm{s}
   = \chi_{ij} w(r_{ij})
   \exp(-\frac{\kappa_j}{\beta} (r_{ij} - \beta s_{0j}))

.. math::
   \dot{\sigma}_{2ij}^\mathrm{o}
   = \chi_{ji} w(r_{ij})
   \exp(-\frac{\kappa_i}{\beta} (r_{ij} - \beta s_{0i}))

and further for unary perfect fcc systems,

Forces
^^^^^^

The forces on atom :math:`i` can be computed as

.. math::
   \mathbf{F}_{i}
   = -\nabla_i E
   = \sum_j \frac{\partial E}{\partial r_{ij}} \frac{\mathbf{r}_{ij}}{r_{ij}}
   = \sum_j \mathbf{f}_{ij}

where the force applied on atom :math:`i` by atom :math:`j` is given by

.. math::
   \mathbf{f}_{ij}
   = \frac{\partial E}{\partial r_{ij}} \frac{\mathbf{r}_{ij}}{r_{ij}}

The derivative of :math:`E` with respect to :math:`r_{ij}` is further written
as

.. math::
   \frac{\partial E}{\partial r_{ij}} = \left(
   \frac{\partial E_{\mathrm{c},i}}{\partial r_{ij}} +
   \frac{\partial E_{\mathrm{c},j}}{\partial r_{ij}} +
   \frac{\partial E_{\mathrm{AS},i}^2}{\partial r_{ij}} +
   \frac{\partial E_{\mathrm{AS},j}^2}{\partial r_{ij}} \right) + \left(
   \frac{\partial E_{\mathrm{AS},i}^1}{\partial r_{ij}} +
   \frac{\partial E_{\mathrm{AS},j}^1}{\partial r_{ij}} \right)

Be careful that we also need to consider the contribution of the energy term
associated to atom :math:`j`.

The first terms depend on :math:`r_{ij}` indirectly via `s_{i}` and `s_{j}`.

.. math::
   \frac{\partial E_{\mathrm{c},i}}{\partial r_{ij}} =
   \frac{\partial E_{\mathrm{c},i}}{\partial s_{i}}
   \frac{\partial s_i}{\partial \sigma_{1i}}
   \frac{\partial \sigma_{1i}}{\partial r_{ij}}

.. math::
   \frac{\partial E_{\mathrm{c},j}}{\partial r_{ij}} =
   \frac{\partial E_{\mathrm{c},j}}{\partial s_{j}}
   \frac{\partial s_j}{\partial \sigma_{1j}}
   \frac{\partial \sigma_{1j}}{\partial r_{ij}}

.. math::
   \frac{\partial E_{\mathrm{AS},i}^2}{\partial r_{ij}} =
   \frac{\partial E_{\mathrm{AS},i}^2}{\partial s_{i}}
   \frac{\partial s_i}{\partial \sigma_{1i}}
   \frac{\partial \sigma_{1i}}{\partial r_{ij}}

.. math::
   \frac{\partial E_{\mathrm{AS},j}^2}{\partial r_{ij}} =
   \frac{\partial E_{\mathrm{AS},j}^2}{\partial s_{j}}
   \frac{\partial s_j}{\partial \sigma_{1j}}
   \frac{\partial \sigma_{1j}}{\partial r_{ij}}

They can be computed using

.. math::
   \frac{\partial E_{\mathrm{c},i}}{\partial s_i}
   = - E_{0i} \lambda_i^2 \dot{s}_i \exp(-\lambda_i \dot{s}_i)

.. math::
   \frac{\partial E_\mathrm{AS}^2}{\partial s_i}
   = -6 V_{0i} \kappa_i \exp(-\kappa_i \dot{s}_i)

.. math::
   \frac{\mathrm{d}s_i}{\mathrm{d}\sigma_{1i}}
   = \frac{\mathrm{d}\dot{s}_i}{\mathrm{d}\sigma_{1i}}
   = -\frac{1}{\beta\eta_{2i}} \frac{1}{\sigma_{1i}}

.. math::
   \frac{\partial \sigma_{1i}}{\partial r_{ij}}
   = \chi_{ij}
   \left(
      \frac{\partial w}{\partial r_{ij}}
      \exp(-\eta_{2j} (r_{ij} - \beta s_{0j})) -
      w(r_{ij}) \eta_{2j}
      \exp(-\eta_{2j} (r_{ij} - \beta s_{0j}))
   \right)
   = \left(
      \frac{1}{w(r_{ij})} \frac{\partial w}{\partial r_{ij}} -
      \eta_{2j}
   \right) \dot{\sigma}_{1ij}^\mathrm{s}

.. math::
   \frac{\partial \sigma_{1j}}{\partial r_{ij}}
   = \chi_{ji}
   \left(
      \frac{\partial w}{\partial r_{ij}}
      \exp(-\eta_{2i} (r_{ij} - \beta s_{0i})) -
      w(r_{ij}) \eta_{2i}
      \exp(-\eta_{2i} (r_{ij} - \beta s_{0i}))
   \right)
   = \left(
      \frac{1}{w(r_{ij})} \frac{\partial w}{\partial r_{ij}} -
      \eta_{2i}
   \right) \dot{\sigma}_{1ij}^\mathrm{o}

The second part directly depends on :math:`r_{ij}` and given by

.. math::
   \frac{\partial E_{\mathrm{AS},i}^1}{\partial r_{ij}} +
   \frac{\partial E_{\mathrm{AS},j}^1}{\partial r_{ij}}
   = - \frac{1}{2} \left(
      \frac{V_{0i}}{2 \gamma_{2i}} \frac{\partial \dot{\sigma}_{2ij}^\mathrm{s}}{\partial r_{ij}} +
      \frac{V_{0j}}{2 \gamma_{2j}} \frac{\partial \dot{\sigma}_{2ij}^\mathrm{o}}{\partial r_{ij}}
   \right)

where

.. math::
   \frac{\partial \dot{\sigma}_{2ij}^\mathrm{s}}{\partial r_{ij}}
   = \chi_{ij}
   \left(
      \frac{\partial w}{\partial r_{ij}}
      \exp(-\frac{\kappa_j}{\beta} (r_{ij} - \beta s_{0j})) -
      w(r_{ij}) \frac{\kappa_j}{\beta}
      \exp(-\frac{\kappa_j}{\beta} (r_{ij} - \beta s_{0j}))
   \right)
   = \left(
      \frac{1}{w(r_{ij})} \frac{\partial w}{\partial r_{ij}} -
      \frac{\kappa_j}{\beta}
   \right)
   \dot{\sigma}_{2ij}^\mathrm{s}

.. math::
   \frac{\partial \dot{\sigma}_{2ij}^\mathrm{o}}{\partial r_{ij}}
   = \chi_{ji}
   \left(
      \frac{\partial w}{\partial r_{ij}}
      \exp(-\frac{\kappa_i}{\beta} (r_{ij} - \beta s_{0i})) -
      w(r_{ij}) \frac{\kappa_i}{\beta}
      \exp(-\frac{\kappa_i}{\beta} (r_{ij} - \beta s_{0i}))
   \right)
   = \left(
      \frac{1}{w(r_{ij})} \frac{\partial w}{\partial r_{ij}} -
      \frac{\kappa_i}{\beta}
   \right)
   \dot{\sigma}_{2ij}^\mathrm{o}

Note that

.. math::
   \frac{\mathrm{d}w}{\mathrm{d}r} = a w(r) (w(r) - 1)

Stress
^^^^^^

The static part of the virial stress can be given as

.. math::
   \tau^{\alpha \beta}
   = \frac{1}{\Omega} \frac{1}{2}
   \sum_{i=1}^{N} \sum_{j \neq i} r_{ij}^{\alpha} f_{ij}^{\beta}
   = \frac{1}{\Omega}
   \sum_{i=1}^{N} \sum_{j > i} r_{ij}^{\alpha} f_{ij}^{\beta}

where :math:`\alpha` and :math:`\beta` are indices for Cartesian components.
When considering all the neighbors for each atom, we should not forget the
factor :math:`1/2`.

Tips
^^^^

For the fcc structure, the numbers of neighbor sites and the
distances of first several shells are

.. math::

   :wowrap:

   n^\mathrm{1NN} &= 12,           & \quad d^\mathrm{1NN} &= \beta s_{i} \\
   n^\mathrm{2NN} &= \phantom{0}6, & \quad d^\mathrm{2NN} &= \sqrt{2}\,d^\mathrm{1NN} \\
   n^\mathrm{3NN} &= 24,           & \quad d^\mathrm{3NN} &= \sqrt{3}\,d^\mathrm{1NN} \\
   n^\mathrm{4NN} &= 12,           & \quad d^\mathrm{4NN} &= \sqrt{4}\,d^\mathrm{1NN} \\
   n^\mathrm{5NN} &= 24,           & \quad d^\mathrm{5NN} &= \sqrt{5}\,d^\mathrm{1NN}

where :math:`s_{i}` is the Wigner–Seitz radius of the species of atom
:math:`i` and :math:`\beta = 2^{-1/2} (16 \pi / 3)^{1/3} \approx 1.809`.
