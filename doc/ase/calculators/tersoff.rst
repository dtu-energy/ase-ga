.. module:: ase.calculators.tersoff
    :synopsis: Tersoff Interatomic Potential

Tersoff Calculator
==================

The Tersoff potential is a three-body empirical potential that can model covalently
bonded solids like silicon or carbon. This implementation provides a native ASE
calculator that follows the `LAMMPS-style Tersoff
<https://docs.lammps.org/pair_tersoff.html>`_ parameterization.

.. note::

    The performance of the routines for this calculator have not been optimized for
    speed nor benchmarked with the LAMMPS implementation.

Theory
------

The many-body interaction is based on bond-order concept that includes both two-body and
three-body terms. The total energy is given by:

.. math::

    \begin{split}
    E & = \frac{1}{2} \sum_i \sum_{j \neq i} V_{ij} \\
    V_{ij} & = f_C(r_{ij}) \left[ f_R(r_{ij}) + b_{ij} f_A(r_{ij}) \right]
    \end{split}

where the repulsive and attractive pair terms are:

.. math::

    \begin{split}
    f_R(r) & = A \exp (-\lambda_1 r) \\
    f_A(r) & = -B \exp (-\lambda_2 r)
    \end{split}

The cutoff function :math:`f_C(r)` ensures smooth decay of interactions:

.. math::

    f_C(r) = \begin{cases}
    1 & r < R - D \\
    \frac{1}{2} - \frac{1}{2} \sin \left( \frac{\pi}{2} \frac{r-R}{D} \right) & R-D < r < R + D \\
    0 & r > R + D
    \end{cases}

The bond order term :math:`b_{ij}` captures the short-range local atomic environment:

.. math::

    \begin{split}
    b_{ij} & = \left( 1 + \beta^n {\zeta_{ij}}^n \right)^{-\frac{1}{2n}} \\
    \zeta_{ij} & = \sum_{k \neq i,j} f_C(r_{ik}) g(\theta_{ijk})
                  \exp \left[ {\lambda_3}^m (r_{ij} - r_{ik})^m \right]
    \end{split}

where :math:`\theta_{ijk}` is the angle between bonds :math:`ij` and :math:`ik`. The
angular function :math:`g(\theta)` is:

.. math::

    g(\theta) = \gamma \left( 1 + \frac{c^2}{d^2} -
                 \frac{c^2}{d^2 + (h - \cos \theta)^2} \right)

where :math:`h = \cos \theta_0` defines the preferred angle :math:`\theta_0`.

The parameters :math:`A`, :math:`B`, :math:`\lambda_1`, :math:`\lambda_2`,
:math:`\lambda_3`, :math:`\beta`, :math:`n`, :math:`c`, :math:`d`, :math:`h`,
:math:`\gamma`, :math:`m`, :math:`R`, and :math:`D` define the potential for each
interaction type.

For a complete description of the functional forms and parameters, see the
:mod:`ase.calculators.tersoff` module documentation.

Parameters
----------

The Tersoff potential is defined by 14 parameters for each three-body interaction:

========= ===================================
Parameter Description
========= ===================================
A         Repulsive pair potential prefactor
B         Attractive pair potential prefactor
lambda1   Decay length for repulsive term
lambda2   Decay length for attractive term
lambda3   Decay length for angular term
beta      Angular strength parameter
n         Angular exponent
c         Angular coefficient
d         Angular parameter
h         Cosine of angle parameter
gamma     Angular scaling
m         Bond order exponent
R         Cutoff distance
D         Cutoff width
========= ===================================

See :class:`ase.calculators.tersoff.TersoffParameters` for details.

Examples
--------

Silicon Diamond Structure
~~~~~~~~~~~~~~~~~~~~~~~~~

Using calculator on silicon crystal in the diamond structure:

.. code-block:: python

    from ase.build import bulk
    from ase.calculators.tersoff import Tersoff, TersoffParameters

    # Create silicon diamond structure
    si = bulk("Si", "diamond", a=5.43)

    # Define Tersoff parameters for Si
    si_params = {
        ("Si", "Si", "Si"): TersoffParameters(
            A=3264.7,
            B=95.373,
            lambda1=3.2394,
            lambda2=1.3258,
            lambda3=1.3258,
            beta=0.33675,
            gamma=1.00,
            m=3.00,
            n=22.956,
            c=4.8381,
            d=2.0417,
            h=0.0000,
            R=3.00,
            D=0.20,
        )
    }

    # Set up calculator
    calc = Tersoff(si_params)
    si.calc = calc

    # Calculate properties
    energy = si.get_potential_energy()
    forces = si.get_forces()
    stress = si.get_stress()

Parameter Updates
~~~~~~~~~~~~~~~~~

The calculator parameters can be updated after initialization:

.. code-block:: python

    # Update single parameters
    calc.set_parameters(("Si", "Si", "Si"), R=2.9, D=0.25)

    # Or replace entire parameter set
    new_params = TersoffParameters(...)
    calc.set_parameters(("Si", "Si", "Si"), params=new_params)

Interface to LAMMPS Files
~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

    This interface has only been tested with `Si.tersoff
    <https://github.com/lammps/lammps/blob/ea1607f1d8a70941cc675f18c4255c25c95c459f/potentials/Si.tersoff>`_
    and `SiC.tersoff
    <https://github.com/lammps/lammps/blob/ea1607f1d8a70941cc675f18c4255c25c95c459f/potentials/SiC.tersoff>`_
    LAMMPS files so it is not necessarily guaranteed to parse all LAMMPS Tersoff files. Therefore, it is
    recommended to format the LAMMPS Tersoff file using similar to the *Si.tersoff* or
    *SiC.tersoff* files.

Read parameters from a `LAMMPS-style Tersoff file
<https://docs.lammps.org/pair_tersoff.html>`_:

.. code-block:: python

    # Initialize from LAMMPS file
    calc = Tersoff.from_lammps("SiC.tersoff")

Tersoff Calculator Class
++++++++++++++++++++++++

.. autoclass:: ase.calculators.tersoff.TersoffParameters
    :class-doc-from: class

.. autoclass:: ase.calculators.tersoff.Tersoff
    :members: from_lammps, read_lammps_format, set_parameters

References
----------

1. J. Tersoff, "New empirical approach for the structure and energy of covalent
      systems", *Phys. Rev. B* **37**, 6991 (1988)
2. J. Tersoff, "Modeling solid-state chemistry: Interatomic potentials for
      multicomponent systems", *Phys. Rev. B* **39**, 5566(R) (1989)
3. S. J. Plimpton, A. Kohlmeyer, A. P. Thompson, et al., "LAMMPS: Large-scale
         Atomic/Molecular Massively Parallel Simulator". Zenodo, Aug. 02, 2023.
         `doi:10.5281/zenodo.10806852 <https://doi.org/10.5281/zenodo.10806852>`_.
