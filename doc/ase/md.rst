==================
Molecular dynamics
==================

.. module:: ase.md
   :synopsis: Molecular Dynamics

Typical computer simulations involve moving the atoms around, either
to optimize a structure (energy minimization) or to do molecular
dynamics.  This chapter discusses molecular dynamics, energy
minimization algorithms will be discussed in the :mod:`ase.optimize`
section.

A molecular dynamics object will operate on the atoms by moving them
according to their forces - it integrates Newton's second law
numerically.  A typical molecular dynamics simulation will use the
`Velocity Verlet dynamics`_.  You create the
:class:`ase.md.verlet.VelocityVerlet` object, giving it the atoms and a time
step, and then you perform dynamics by calling its
:meth:`~verlet.VelocityVerlet.run` method::

  dyn = VelocityVerlet(atoms, dt=5.0 * units.fs,
                       trajectory='md.traj', logfile='md.log')
  dyn.run(1000)  # take 1000 steps

A number of algorithms can be used to perform molecular
dynamics, with slightly different results.


.. note::

   Prior to ASE version 3.21.0, inconsistent units were used to
   specify temperature.  Some modules expected kT (in eV), others T
   (in Kelvin).  From ASE 3.21.0, all molecular dynamics modules
   expecting a temperature take a parameter ``temperature_K`` which is
   the temperature in Kelvin.  For compatibility, they still accept
   the ``temperature`` parameter in the same unit as previous versions
   (eV or K), but using the old parameter will issue a warning.



Choosing the time step
======================

All the dynamics objects need a time step.  Choosing it too small will
waste computer time, choosing it too large will make the dynamics
unstable, typically the energy increases dramatically (the system
"blows up").  If the time step is only a little to large, the lack of
energy conservation is most obvious in `Velocity Verlet dynamics`_,
where energy should otherwise be conserved.

Experience has shown that 5 femtoseconds is a good choice for most metallic
systems.  Systems with light atoms (e.g. hydrogen) and/or with strong
bonds (carbon) will need a smaller time step.

All the dynamics objects documented here are sufficiently related to
have the same optimal time step.


File output
===========

The time evolution of the system can be saved in a trajectory file,
by creating a trajectory object, and attaching it to the dynamics
object.  This is documented in the module :mod:`ase.io.trajectory`.
You can attach the trajectory explicitly to the dynamics object, and
you may want to use the optional ``interval`` argument, so every
time step is not written to the file.

Alternatively, you can just use the ``trajectory`` keyword when
instantiating the dynamics object as in the example above. In this
case, a ``loginterval`` keyword may also be supplied to specify the
frequency of writing to the trajectory. The loginterval keyword will
apply to both the trajectory and the logfile.


Logging
=======

A logging mechanism is provided, printing time; total, potential and
kinetic energy; and temperature (calculated from the kinetic energy).
It is enabled by giving the ``logfile`` argument when the dynamics
object is created, ``logfile`` may be an open file, a filename or the
string '-' meaning standard output.  Per default, a line is printed
for each timestep, specifying the ``loginterval`` argument will chance
this to a more reasonable frequency.

The logging can be customized by explicitly attaching a
:class:`MDLogger` object to the dynamics::

  from ase.md import MDLogger
  dyn = VelocityVerlet(atoms, dt=2*ase.units.fs)
  dyn.attach(MDLogger(dyn, atoms, 'md.log', header=False, stress=False,
             peratom=True, mode="a"), interval=1000)

This example will skip the header line and write energies per atom
instead of total energies.  The parameters are

  ``header``: Print a header line defining the columns.

  ``stress``: Print the six components of the stress tensor.

  ``peratom``:  Print energy per atom instead of total energy.

  ``mode``:  If 'a', append to existing file, if 'w' overwrite
  existing file.

Despite appearances, attaching a logger like this does *not* create a
cyclic reference to the dynamics.

.. note::

   If building your own logging class, be sure not to attach the dynamics
   object directly to the logging object. Instead, create a weak reference
   using the ``proxy`` method of the ``weakref`` package. See the
   *ase.md.MDLogger* source code for an example. (If this is not done, a
   cyclic reference may be created which can cause certain calculators
   to not terminate correctly.)


.. autoclass:: MDLogger


Constant NVE simulations (the microcanonical ensemble)
======================================================

Newton's second law preserves the total energy of the system, and a
straightforward integration of Newton's second law therefore leads to
simulations preserving the total energy of the system (E), the number
of atoms (N) and the volume of the system (V).  The most appropriate
algorithm for doing this is velocity Verlet dynamics, since it gives
very good long-term stability of the total energy even with quite
large time steps.  Fancier algorithms such as Runge-Kutta may give
very good short-term energy preservation, but at the price of a slow
drift in energy over longer timescales, causing trouble for long
simulations.

In a typical NVE simulation, the temperature will remain approximately
constant, but if significant structural changes occurs they may result
in temperature changes.  If external work is done on the system, the
temperature is likely to rise significantly.


Velocity Verlet dynamics
------------------------

.. module:: ase.md.verlet

.. autoclass:: VelocityVerlet


``VelocityVerlet`` is the only dynamics implementing the NVE ensemble.
It requires two arguments, the atoms and the time step.  Choosing
a too large time step will immediately be obvious, as the energy will
increase with time, often very rapidly.

Example: See the tutorial :ref:`md_tutorial`.


Constant NVT simulations (the canonical ensemble)
=================================================

Since Newton's second law conserves energy and not temperature,
simulations at constant temperature will somehow involve coupling the
system to a heat bath.  This cannot help being somewhat artificial.
Two different approaches are possible within ASE.  In Langevin
dynamics, each atom is coupled to a heat bath through a fluctuating
force and a friction term.  In Nosé-Hoover dynamics, a term
representing the heat bath through a single degree of freedom is
introduced into the Hamiltonian.


Langevin dynamics
-----------------

.. module:: ase.md.langevin

.. autoclass:: Langevin

The Langevin class implements Langevin dynamics, where a (small)
friction term and a fluctuating force are added to Newton's second law
which is then integrated numerically.  The temperature of the heat
bath and magnitude of the friction is specified by the user, the
amplitude of the fluctuating force is then calculated to give that
temperature.  This procedure has some physical justification: in a
real metal the atoms are (weakly) coupled to the electron gas, and the
electron gas therefore acts like a heat bath for the atoms.  If heat
is produced locally, the atoms locally get a temperature that is
higher than the temperature of the electrons, heat is transferred to
the electrons and then rapidly transported away by them.  A Langevin
equation is probably a reasonable model for this process.

A disadvantage of using Langevin dynamics is that if significant heat
is produced in the simulation, then the temperature will stabilize at
a value higher than the specified temperature of the heat bath, since
a temperature difference between the system and the heat bath is
necessary to get a finite heat flow.  Another disadvantage is that the
fluctuating force is stochastic in nature, so repeating the simulation
will not give exactly the same trajectory.

When the ``Langevin`` object is created, you must specify a time step,
a temperature (in Kelvin) and a friction coefficient.
A typical range for the friction coefficient may be 0.001–0.1 fs\ :sup:`−1`
(1–100 ps\ :sup:`−1`).
For example, you can give a friction coefficient of 0.01 fs\ :sup:`−1`
(10 ps\ :sup:`−1`) as

.. code-block:: python

  dyn = Langevin(
      atoms,
      timestep=5.0 * units.fs,
      temperature_K=300.0,  # temperature in K
      friction=0.01 / units.fs,
  )

Both the friction and the temperature can be replaced with arrays
giving per-atom values.  This is mostly useful for the friction, where
one can choose a rather high friction near the boundaries, and set it
to zero in the part of the system where the phenomenon being studied
is located.


Andersen dynamics
-----------------

.. module:: ase.md.andersen

.. autoclass:: Andersen

The Andersen class implements Andersen dynamics, where constant
temperature is imposed by stochastic collisions with a heat bath.
With a (small) probability (`andersen_prob`) the collisions act
occasionally on velocity components of randomly selected particles
Upon a collision the new velocity is drawn from the
Maxwell-Boltzmann distribution at the corresponding temperature.
The system is then integrated numerically at constant energy
according to the Newtonian laws of motion. The collision probability
is defined as the average number of collisions per atom and timestep. 
The algorithm generates a canonical distribution. [1] However, due
to the random decorrelation of velocities, the dynamics are
unphysical and cannot represent dynamical properties like e.g.
diffusion or viscosity. Another disadvantage is that the collisions
are stochastic in nature, so repeating the simulation will not give
exactly the same trajectory.

When the ``Andersen`` object is created, you must specify a time step,
a temperature (in Kelvin) and a collision probability. Typical 
values for this probability are in the order of 1e-4 to 1e-1.

::

  # Room temperature simulation (300 Kelvin, Andersen probability: 0.002)
  dyn = Andersen(atoms, 5 * units.fs, 300, 0.002)

References:

[1] D. Frenkel and B. Smit, Understanding Molecular Simulation
(Academic Press, London, 1996)


Nosé-Hoover dynamics
--------------------

In Nosé-Hoover dynamics, an extra term is added to the Hamiltonian
representing the coupling to the heat bath.  From a pragmatic point of
view one can regard Nosé-Hoover dynamics as adding a friction term to
Newton's second law, but dynamically changing the friction coefficient
to move the system towards the desired temperature.  Typically the
"friction coefficient" will fluctuate around zero.

Nosé-Hoover dynamics is not implemented as a separate class, but is a
special case of NPT dynamics.


Berendsen NVT dynamics
-----------------------
.. module:: ase.md.nvtberendsen

.. autoclass:: NVTBerendsen

In Berendsen NVT simulations the velocities are scaled to achieve the desired
temperature. The speed of the scaling is determined by the parameter taut.

This method does not result proper NVT sampling but it usually is
sufficiently good in practice (with large taut). For discussion see
the gromacs manual at www.gromacs.org.

::

  # Room temperature simulation (300K, 0.1 fs time step)
  dyn = NVTBerendsen(atoms, 0.1 * units.fs, 300, taut=0.5*1000*units.fs)



Constant NPT simulations (the isothermal-isobaric ensemble)
===========================================================

.. module:: ase.md.npt

.. autoclass:: NPT

    .. automethod:: run 
    .. automethod:: set_stress
    .. automethod:: set_temperature 
    .. automethod:: set_mask 
    .. automethod:: set_fraction_traceless 
    .. automethod:: get_strain_rate 
    .. automethod:: set_strain_rate 
    .. automethod:: get_time 
    .. automethod:: initialize
    .. automethod:: get_gibbs_free_energy 
    .. automethod:: zero_center_of_mass_momentum 
    .. automethod:: attach 

       
Berendsen NPT dynamics
-----------------------
.. module:: ase.md.nptberendsen

.. autoclass:: NPTBerendsen

In Berendsen NPT simulations the velocities are scaled to achieve the desired
temperature. The speed of the scaling is determined by the parameter taut.

The atom positions and the simulation cell are scaled in order to achieve
the desired pressure.

This method does not result proper NPT sampling but it usually is
sufficiently good in practice (with large taut and taup). For discussion see
the gromacs manual at www.gromacs.org. or amber at ambermd.org


::

  # Room temperature simulation (300K, 0.1 fs time step, atmospheric pressure)
  dyn = NPTBerendsen(atoms, timestep=0.1 * units.fs, temperature_K=300,
                     taut=100 * units.fs, pressure_au=1.01325 * units.bar,
                     taup=1000 * units.fs, compressibility_au=4.57e-5 / units.bar)

Bussi dynamics
--------------

.. module:: ase.md.bussi

.. autoclass:: Bussi

The Bussi class implements the Bussi dynamics, where constant
temperature is imposed by a stochastic velocity rescaling algorithm.
The thermostat is conceptually very close to the Berendsen thermostat,
but does sample the canonical ensemble correctly. Given that the thermostat
is both correct and global, it is advised to use it for production runs.

Contour Exploration
-------------------
.. module:: ase.md.contour_exploration

.. autoclass:: ContourExploration

Contour Exploration evolves the system along constant potentials energy
contours on the potential energy surface. The method uses curvature based
extrapolation and a potentiostat to correct for potential energy errors. It is
similar to molecular dynamics but with a potentiostat rather than a thermostat.
Without changes in kinetic energy, it is more useful to automatically scale
step sizes to the curvature of the potential energy contour via an
``angle_limit`` while enforcing a ``maxstep`` to ensure potentiostatic 
accuracy. [1] To escape loops on the pontential energy surface or to break
symmetries, a random drift vector parallel to the contour can be applied as a
fraction of the step size via ``parallel_drift``. Contour exploration cannot 
be used at minima since the contour is a single point.


::

  # Contour exploration at the current potential energy
  dyn = ContourExploration(atoms)

References:

[1] M. J. Waters and J. M. Rondinelli, `Contour Exploration with Potentiostatic
Kinematics` ArXiv:2103.08054 (https://arxiv.org/abs/2103.08054)



Velocity distributions
======================

A selection of functions are provided to initialize atomic velocities
to the correct temperature.

.. module:: ase.md.velocitydistribution

.. autofunction:: MaxwellBoltzmannDistribution

.. autofunction:: Stationary

.. autofunction:: ZeroRotation

.. autofunction:: PhononHarmonics

.. autofunction:: phonon_harmonics

Post-simulation Analysis
========================

Functionality is provided to perform analysis of atomic/molecular behaviour as calculation in a molecular dynamics simulation. Currently, this is presented as a class to address the Einstein equation for diffusivity.

.. module:: ase.md.analysis

.. autoclass:: DiffusionCoefficient

