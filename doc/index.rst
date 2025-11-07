.. ASE‑GA documentation master file

=========================
ASE‑GA: Genetic Algorithm
=========================

Welcome to the **ASE‑GA** documentation!

ASE‑GA (Atomic Simulation Environment – Genetic Algorithm) was originally implemented within the ASE package for global optimization of atomic structures using genetic algorithms. It has since been spun off into a standalone package, enabling more modular development and focused growth while also enabling ASE to focus on core capabilities.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Tutorials

   tutorials/ga_optimize.rst
   tutorials/ga_fcc_alloys.rst
   tutorials/ga_convex_hull.rst
   tutorials/ga_bulk.rst
   tutorials/ga_molecular_crystal.rst
   api.rst

Overview
--------

Genetic algorithms (GAs) offer an effective alternative to Monte Carlo-style search methods for finding global minima. ASE‑GA implements a GA framework allowing users to define populations, pairing operators, mutations, and stopping conditions tailored to their optimization problem within materials and atomic configurations. The method is described in detail in the following publications, that also shows the large variability of systems it is able to handle:

For **small clusters on/in support material** in:

   | L. B. Vilhelmsen and B. Hammer
   | :doi:`A genetic algorithm for first principles global structure optimization of supported nano structures <10.1063/1.4886337>`
   | The Journal of chemical physics, Vol. 141 (2014), 044711

For **medium sized alloy clusters** in:

   | S. Lysgaard, D. D. Landis, T. Bligaard and T. Vegge
   | :doi:`Genetic Algorithm Procreation Operators for Alloy Nanoparticle Catalysts <10.1007/s11244-013-0160-9>`
   | Topics in Catalysis, Vol **57**, No. 1-4, pp. 33-39, (2014)
   
A search for **mixed metal ammines for ammonia storage** have been performed
using the GA in:

   | P. B. Jensen, S. Lysgaard, U. J. Quaade and T. Vegge
   | :doi:`Designing Mixed Metal Halide Ammines for Ammonia Storage Using Density Functional Theory and Genetic Algorithms <10.1039/C4CP03133D>`
   | Physical Chemistry Chemical Physics, Vol **16**, No. 36, pp. 19732-19740, (2014)

**Crystal structures searches** was performed in:

   | M. Van den Bossche, H. Grönbeck, and B. Hammer
   | :doi:`Tight-Binding Approximation-Enhanced Global Optimization <10.1021/acs.jctc.8b00039>`
   | J. Chem. Theory Comput. 2018, **14**, 2797−2807
   
   
Installation
------------

Install the **ASE‑GA** package via pip or from source:

.. code-block:: bash

   pip install ase-ga

Or clone and install from GitHub:

.. code-block:: bash

   git clone https://github.com/dtu-energy/ase-ga.git
   cd ase-ga
   pip install -e .

..
   Quickstart
   ----------

   Here's a minimal example to initialize and run a genetic algorithm search:

   .. code-block:: python

      from ase_ga.data import PrepareDB

      # Setup database and initial population
      PrepareDB(...)
      StartGenerator(...)

      # Run optimizer

   See the detailed tutorial ``tutorial_ga.rst`` for full examples on cluster search, bulk crystal optimization, and running in parallel.

Tutorials
---------

- :ref:`genetic_algorithm_optimization_tutorial`: Covers the different aspects of running the GA both locally and on an HPC with a cluster on a surface example.
- :ref:`fcc_alloys_tutorial`: Outlines a GA search with an example of finding the most stable fcc alloys.
- :ref:`convex_hull_tutorial`: Determine the full convex hull in a single GA run. It introduces ranked population that group candidates according to a variable.
- :ref:`ga_bulk_tutorial`: Describes a crystal structure search and introduces bulk-specific operators like `StrainMutation`, `PermuStrainMutation`, and `SoftMutation`
- :ref:`ga_molecular_crystal_tutorial`: Search for molecular crystals. It resembles the crystal structure search but keeps molecular identity of the molecules.

Implementation and API Reference
--------------------------------

The GA implementation is diverse. It is structured such that it can be tailored to the specific problem investigated and to the computational resources available (single computer or a large computer cluster).

A genetic algorithm (GA) applies a Darwinian approach to optimization by maintaining a :doc:`Population <api/ase_ga.population>` of solution `Candidates` (always represented with :py:class:`Atoms <ase:ase.Atoms>`) to a given problem. Over successive iterations, the population evolves toward better solutions through :doc:`Mating and Mutation <api/ase_ga.offspring_creator>` of selected candidates, with the fittest resulting offspring added back into the population (as new candidates).

The :doc:`raw score <api/ase_ga.raw_score>` of a candidate is determined by a function that evaluates its quality - for example, by measuring the stability or performance of a structure. The fitness is calculated when the candidate enters the population, it is a measure of how well the candidate performs compared with the rest of the population. To keep the population size constant, the least fit candidates are removed.

Mating (or :ref:`Crossover <crossover_api>`) combines multiple candidates to produce new offspring that inherit traits from their parents. When beneficial traits are passed on and combined, the population improves. However, relying only on crossover can reduce diversity: repeatedly mating very similar candidates tends to stagnate progress.  

:class:`Mutation` introduces diversity by making random alterations to candidates, preventing premature convergence and helping the algorithm explore a wider search space. :class:`Comparators` prevent identical (or just too similar candidates) in being present in the population, thereby also enhancing diversity.

A GA run is at an end when no evolution is identified in the population for a certain period. See the :doc:`Convergence` for the available parameters.

See the :ref:`api` for the full auto-generated API docs.

Examples
--------

Coming up

Changelog
---------

Notably, ASE version 3.9.0 (May 2015) first introduced the ``ase.ga`` module, and recent versions have extended support for crystal structure prediction operators.

Releases include:

- **v1.0.0** – Initial extraction from ASE  

Contributing
------------

Contributions are welcome! Please follow the coding style, tests, and guidelines in :doc:`contributing`.

Support
-------

If you have questions or encounter issues:

- File an issue at the project’s GitHub repository https://github.com/dtu-energy/ase-ga/issues

