ASE ecosystem
=============

This is a list of software packages related to ASE or using ASE.
These could well be of interest to ASE users in general.
If you know of a project which
should be listed here, but isn't, please open a merge request adding
link and descriptive paragraph.

Listed in alphabetical order, for want of a better approach.

 * `abTEM <https://abtem.readthedocs.io/en/latest/index.html>`_:
   abTEM provides a Python API for running simulations of (scanning)
   transmission electron microscopy images and diffraction patterns.

 * `ACAT <https://asm-dtu.gitlab.io/acat/>`_:
   ACAT is a Python package for atomistic modelling of metal or alloy 
   heterogeneoues catalysts. ACAT provides automatic identification of 
   adsorption sites and adsorbate coverages for a wide range of surfaces 
   and nanoparticles. ACAT also provides tools for structure generation 
   and global optimization of catalysts with and without adsorbates.

 * `AGOX <https://gitlab.com/agox/agox/>`_:
   The Atomistic Global Optimization X package contains a collection of 
   tools for global optimization of atomic systems. The package allows 
   running a variety of standard global optimization algorithms, such as random structure
   search, basin-hopping in addition to machine-learning enhanced algorithms like 
   GOFEE. Any ASE calculator can be used as the objective function for the optimization.

 * `atomicrex <https://atomicrex.org/>`_:
   atomicrex is a versatile tool for the construction of interatomic
   potential models. It includes a Python interface for integration
   with first-principles codes via ASE as well as other Python
   libraries.

 * `CHGNet <https://github.com/CederGroupHub/chgnet>`_:
   A pretrained universal neural network potential for charge-informed
   atomistic modeling

 * `CLEASE <https://gitlab.com/computationalmaterials/clease#clease>`_:
   CLuster Expansion in Atomic Simulation Environment (CLEASE) is a package
   that automates the cumbersome setup and construction procedure of cluster
   expansion (CE). It provides a comprehensive list of tools for specifying
   parameters for CE, generating training structures, fitting effective cluster
   interaction (ECI) values and running Monte Carlo simulations.

 * `COGEF <https://cogef.gitlab.io/cogef/>`_:
   COnstrained Geometries simulate External Force.  This
   package is useful for analysing properties of bond-breaking
   reactions, such as how much force is required to break a chemical
   bond.

 * `DebyeCalculator <https://github.com/FrederikLizakJohansen/DebyeCalculator>`_:
   A vectorised implementation of the Debye Scattering Equation on CPU and GPU to calculate the scattering intensity I(Q), the Total Scattering Structure
   Function S(Q), the Reduced Total Scattering Function F(Q), or the Reduced Atomic Pair Distribution Function G(r) from an atomic structure. Use 
   DebyeCalculator to simulate powder diffraction, total scattering with pair distribution function or small-angle scattering data of finite systems such as   
   nanoparticles.

 * `effmass <https://github.com/lucydot/effmass/>`_:
   Calculates various definitions of effective mass from the electronic 
   bandstructure of a semiconductor.

 * `evgraf <https://github.com/pmla/evgraf>`_:
   A python library for crystal reduction (i.e. finding primitive cells), and
   identification and symmetrization of structures with inversion
   pseudosymmetry.

 * `FHI-vibes <https://vibes-developers.gitlab.io/vibes/>`_:
   A python package for calculating and analyzing the vibrational properties
   of solids from first principles. FHI-vibes bridges between the harmonic
   approximation and fully anharmonic molecular dynamics simulations.
   FHI-vibes builds on several existing packages including ASE, and provides
   a consistent and user-friendly interface.

 * `gpatom <https://gitlab.com/gpatom/ase-gpatom>`_: APython package
   which provides several tools for geometry optimisation and related
   tasks in atomistic systems using machine learning surrogate models.
   gpatom is an extension to the Atomic Simulation Environment.

 * `hiphive <https://hiphive.materialsmodeling.org>`_:
   hiPhive is a tool for efficiently extracting high-order force
   constants. It is interfaced with ASE to enable easy integration
   with first-principles codes. hiphive also provides an ASE-style
   calculator to enable sampling of force constant expansions via
   molecular dynamics simulations.

 * `icet <https://icet.materialsmodeling.org/>`_:
   The integration cluster expansion toolkit. icet is a flexible and
   extendable software package for constructing and sampling alloy
   cluster expansions. It supports a wide range of regression and
   validation techniques, and includes a Monte Carlo module with
   support for many different thermodynamic ensembles.

 * `matgl <https://github.com/materialsvirtuallab/matgl>`_:
   Graph deep learning library for materials

 * `matscipy <https://github.com/libAtoms/matscipy>`_:
   matscipy is a generic materials science toolbox built around ASE.
   It provides useful routines for plasticity and dislocations, fracture
   mechanics, electro-chemistry, tribology, and elastic properties.
   In addition to domain-specific routines, it also implements a set of
   general-purpose, low-level utilities such as efficient neighbour lists.

 * `NequIP <https://github.com/mir-group/nequip>`_:
   Euclidian Equivariant neural network potentials.  Nequip can fit
   neural network potentials to series of DFT calculations (using
   e.g. ASE trajectory files), and then be used to perform
   optimization and molecular dynamics in ASE or LAMMPS.

 * `QuAcc <https://github.com/Quantum-Accelerators/quacc>`_:
   A flexible platform for high-throughput, database-driven computational 
   materials science and quantum chemistry workflows built around ASE.

 * `SchNet Pack <https://github.com/atomistic-machine-learning/schnetpack>`_:
   Deep Neural Networks for Atomistic Systems

 * `Sella <https://github.com/zadorlab/sella>`_:
   Sella is a saddle point refinement (optimization) tool which uses
   the `Optimize <ase/optimize.html>`_ API. Sella supports minimization and
   refinement of arbitrary-order saddle points with constraints.
   Additionally, Sella can perform intrinsic reaction coordinate (IRC)
   calculations.

 * `TorchANI <https://github.com/aiqm/torchani>`_:
   Accurate Neural Network Potential on PyTorch

 * `Wulffpack <https://wulffpack.materialsmodeling.org/>`_:
   Python package for making Wulff constructions, typically for finding
   equilibrium shapes of nanoparticles. WulffPack constructs both continuum
   models and atomistic structures for further modeling with, e.g., molecular
   dynamics or density functional theory.
