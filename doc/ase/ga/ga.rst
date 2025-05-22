.. module:: ase.ga
   :synopsis: Genetic Algorithm Optimization

===================
 Genetic Algorithm
===================

Genetic algorithms (GA) have proven a good alternative to Monte Carlo
type optimization methods for global structure and materials properties optimization. A GA has
recently been implemented into ase.

The use of the GA is best learned through :ref:`tutorials`.

The GA implementation is diverse. It (or previous versions of it) has been used in publications with differing subjects such as structure of gold clusters on surfaces, composition of alloy nanoparticles, ammonia storage in mixed metal ammines and more. The implementation is structured such that it can be tailored to the specific problem investigated and to the computational resources available (single computer or a large computer cluster).
   
The method is described in detail in the following publications:

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
   
A simple tutorial explaining how to set up a database and perform a
similar search can be found here: :ref:`fcc_alloys_tutorial`
