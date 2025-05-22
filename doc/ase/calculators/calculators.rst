.. module:: ase.calculators
   :synopsis: Energy, force and stress calculators.

.. _calculators:

===========
Calculators
===========

For ASE, a calculator is a black box that can take atomic numbers and
atomic positions from an :class:`~ase.Atoms` object and calculate the
energy and forces and sometimes also stresses.

In order to calculate forces and energies, you need to attach a
calculator object to your atoms object:

>>> atoms = read('molecule.xyz')
>>> e = atoms.get_potential_energy()  # doctest: IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/jjmo/ase/atoms/ase.py", line 399, in get_potential_energy
    raise RuntimeError('Atoms object has no calculator.')
RuntimeError: Atoms object has no calculator.
>>> from ase.calculators.abinit import Abinit
>>> calc = Abinit(...)
>>> atoms.calc = calc
>>> e = atoms.get_potential_energy()
>>> print(e)
-42.0

Here we attached
an instance of the :mod:`ase.calculators.abinit` class and then
we asked for the energy.


.. _supported calculators:

Supported calculators
=====================

The calculators can be divided in four groups:

1) Abacus_, ALIGNN_, AMS_, Asap_, BigDFT_, CHGNet_, DeePMD-kit_, DFTD3_, DFTD4_, DFTK_, EquiFormerV2_, FLEUR_, GPAW_, Hotbit_, M3GNet_, MACE_, OrbModels_, SevenNet_, TBLite_, and XTB_
   have their own native or external ASE interfaces.

2) ABINIT, AMBER, CP2K, CASTEP, deMon2k, DFTB+, ELK, EXCITING, FHI-aims, GAUSSIAN,
   Gromacs, LAMMPS, MOPAC, NWChem, Octopus, ONETEP, PLUMED, psi4, Q-Chem, Quantum ESPRESSO, SIESTA,
   TURBOMOLE and VASP, have Python wrappers in the ASE package, but the actual
   FORTRAN/C/C++ codes are not part of ASE.

3) Pure python implementations included in the ASE package: EMT, EAM,
   Lennard-Jones, Morse, Tersoff, and HarmonicCalculator.

4) Calculators that wrap others, included in the ASE package:

   - :class:`ase.calculators.checkpoint.CheckpointCalculator`
   - :class:`ase.calculators.fd.FiniteDifferenceCalculator`
   - :class:`ase.calculators.loggingcalc.LoggingCalculator`
   - :class:`ase.calculators.mixing.LinearCombinationCalculator`
   - :class:`ase.calculators.mixing.MixedCalculator`
   - :class:`ase.calculators.mixing.SumCalculator`
   - :class:`ase.calculators.mixing.AverageCalculator`
   - :class:`ase.calculators.socketio.SocketIOCalculator`

   - :ref:`Grimme-D3 <grimme>` potential

   - QM/MM calculators

     - :class:`~ase.calculators.qmmm.EIQMMM`
     - :class:`~ase.calculators.qmmm.SimpleQMMM`

========================================= ===========================================
name                                      description
========================================= ===========================================
Abacus_                                   DFT supporting both pw and lcao basis
ALIGNN_                                   Atomistic Line Graph Neural Network force field
AMS_                                      Amsterdam Modeling Suite
Asap_                                     Highly efficient EMT code
BigDFT_                                   Wavelet based code for DFT
CHGNet_                                   Universal neural network  potential for charge-informed atomistics
DeePMD-kit_                               A deep learning package for many-body potential energy representation
DFTD3_                                    London-dispersion correction
DFTD4_                                    Charge-dependent London-dispersion correction
DFTK_                                     Plane-wave code for DFT and related models
EquiFormerV2_                             Equivariant graph-based denoising transformer universal potential
FLEUR_                                    Full Potential LAPW code
GPAW_                                     Real-space/plane-wave/LCAO PAW code
Hotbit_                                   DFT based tight binding
M3GNet_                                   Materials 3-body Graph Network universal potential
MACE_                                     Many-body potential using higher-order equivariant message passing
OrbModels_                                Fast, scalable, universal GNN potentials with diffusion pretraining
SevenNet_                                 Scalable EquiVariance Enabled Neural Network interatomic potential
TBLite_                                   Light-weight tight-binding framework
XTB_                                      Semiemprical extended tight-binding program package
:mod:`~ase.calculators.abinit`            Plane-wave pseudopotential code
:mod:`~ase.calculators.amber`             Classical molecular dynamics code
:mod:`~ase.calculators.castep`            Plane-wave pseudopotential code
:mod:`~ase.calculators.cp2k`              DFT and classical potentials
:mod:`~ase.calculators.demon`             Gaussian based DFT code
:mod:`~ase.calculators.demonnano`         DFT based tight binding code
:mod:`~ase.calculators.dftb`              DFT based tight binding
:mod:`~ase.calculators.dmol`              Atomic orbital DFT code
:mod:`~ase.calculators.eam`               Embedded Atom Method
elk                                       Full Potential LAPW code
:mod:`~ase.calculators.espresso`          Plane-wave pseudopotential code
:mod:`~ase.calculators.exciting`          Full Potential LAPW code
:mod:`~ase.calculators.aims`              Numeric atomic orbital, full potential code
:mod:`~ase.calculators.gamess_us`         Gaussian based electronic structure code
:mod:`~ase.calculators.gaussian`          Gaussian based electronic structure code
:mod:`~ase.calculators.gromacs`           Classical molecular dynamics code
:mod:`~ase.calculators.gulp`              Interatomic potential code
:mod:`~ase.calculators.harmonic`          Hessian based harmonic force-field code
:mod:`~ase.calculators.kim`               Classical MD with standardized models
:mod:`~ase.calculators.lammps`            Classical molecular dynamics code
:mod:`~ase.calculators.mixing`            Combination of multiple calculators
:mod:`~ase.calculators.mopac`             Semiempirical molecular orbital code
:mod:`~ase.calculators.nwchem`            Gaussian based electronic structure code
:mod:`~ase.calculators.octopus`           Real-space pseudopotential code
:mod:`~ase.calculators.onetep`            Linear-scaling pseudopotential code
:mod:`~ase.calculators.openmx`            LCAO pseudopotential code
:mod:`~ase.calculators.orca`              Gaussian based electronic structure code
:mod:`~ase.calculators.plumed`            Enhanced sampling method library
:mod:`~ase.calculators.psi4`              Gaussian based electronic structure code
:mod:`~ase.calculators.qchem`             Gaussian based electronic structure code
:mod:`~ase.calculators.siesta`            LCAO pseudopotential code
:mod:`~ase.calculators.turbomole`         Fast atom orbital code
:mod:`~ase.calculators.tersoff`           Tersoff bond-order potential
:mod:`~ase.calculators.vasp`              Plane-wave PAW code
:mod:`~ase.calculators.emt`               Effective Medium Theory calculator
lj                                        Lennard-Jones potential
morse                                     Morse potential
:mod:`~ase.calculators.checkpoint`        Checkpoint calculator
:mod:`~ase.calculators.fd`                Finite-difference calculator
:mod:`~ase.calculators.loggingcalc`       Logging calculator
:mod:`~ase.calculators.socketio`          Socket-based interface to calculators
:mod:`~ase.calculators.dftd3`             DFT-D3 dispersion correction calculator
:class:`~ase.calculators.qmmm.EIQMMM`     Explicit Interaction QM/MM
:class:`~ase.calculators.qmmm.SimpleQMMM` Subtractive (ONIOM style) QM/MM
========================================= ===========================================

.. index:: D3, Grimme
.. _grimme:

.. note::

    A Fortran implemetation of the Grimme-D3 potential, that can be used as
    an add-on to any ASE calculator, can be found here:
    https://gitlab.com/ehermes/ased3/tree/master.

The calculators included in ASE are used like this:

>>> from ase.calculators.abc import ABC
>>> calc = ABC(...)

where ``abc`` is the module name and ``ABC`` is the class name.

.. _Abacus: https://gitlab.com/1041176461/ase-abacus
.. _ALIGNN: https://github.com/usnistgov/alignn?tab=readme-ov-file#alignnff
.. _AMS: https://www.scm.com/doc/plams/examples/AMSCalculator/ASECalculator.html#asecalculatorexample
.. _Asap: https://wiki.fysik.dtu.dk/asap
.. _BigDFT: https://l_sim.gitlab.io/bigdft-suite/tutorials/Interoperability-Simulation.html#ASE-Interoperability
.. _CHGNet: https://github.com/CederGroupHub/chgnet/blob/e2a2b82bf2c64e5a3d39cd75d0addfa864a2771a/chgnet/model/dynamics.py#L63
.. _GPAW: https://gpaw.readthedocs.io
.. _Hotbit: https://github.com/pekkosk/hotbit
.. _DFTK: https://dftk.org
.. _DeePMD-kit: https://github.com/deepmodeling/deepmd-kit
.. _DFTD4: https://github.com/dftd4/dftd4/tree/main/python
.. _DFTD3: https://dftd3.readthedocs.io/en/latest/api/python.html#module-dftd3.ase
.. _EquiFormerV2: https://github.com/FAIR-Chem/fairchem#quick-start
.. _FLEUR: https://github.com/JuDFTteam/ase-fleur
.. _M3GNet: https://matgl.ai/matgl.ext.html#class-matglextasem3gnetcalculatorpotential-potential-state_attr-torchtensor--none--none-stress_weight-float--10-kwargs
.. _MACE: https://mace-docs.readthedocs.io/en/latest/guide/ase.html
.. _OrbModels: https://github.com/orbital-materials/orb-models/tree/main#usage-with-ase-calculator
.. _SevenNet: https://github.com/MDIL-SNU/SevenNet#ase-calculator
.. _TBLite: https://tblite.readthedocs.io/en/latest/users/ase.html
.. _XTB: https://xtb-python.readthedocs.io/en/latest/ase-calculator.html

Calculator configuration
========================

Calculators that depend on external codes or files are generally
configurable.  ASE loads the configuration from a configfile located
at ``~/.config/ase/config.ini``.  The default path can be overriden by
setting the environment variable ``ASE_CONFIG_PATH`` to another path
or paths separated by colon.

To see the full configuration on a given machine, run
:command:`ase info --calculators`.

An example of a config file is as follows::

    [abinit]
    command = mpiexec /usr/bin/abinit
    pp_paths = /usr/share/abinit/pseudopotentials

    [espresso]
    command = mpiexec pw.x
    pseudo_path = /home/ase/upf_pseudos

Calculators build a full command by appending command-line arguments
to the configured command.  Therefore, the command should normally consist
of any parallel arguments followed by the binary, but should not
include further flags unless desired for a specific reason.
The command is also used to build a full command for e.g.
socket I/O calculators.

The Espresso calculator can then invoked in the following way::

    >>> from ase.build import bulk
    >>> from ase.calculators.espresso import Espresso
    >>> espresso = Espresso(
                       input_data = {
                            'system': {
                               'ecutwfc': 60,
                            }},
                       pseudopotentials = {'Si': 'si_lda_v1.uspp.F.UPF'},
                       )
    >>> si = bulk('Si')
    >>> si.calc = espresso
    >>> si.get_potential_energy()
    -244.76638508140397

It can be useful for software libraries to override the local
configuration.  To do so, the code should supply the configurable
information by instantiating a “profile”, e.g.,
``Abinit(profile=AbinitProfile(command=command))``.  The profile
encloses the configurable information specific to a particular code,
so this may differ depending on which code.  It can also be
useful for software libraries that manage their own configuration
to set the ``ASE_CONFIG_PATH`` to an empty string.


.. toctree::

   eam
   emt
   abinit
   amber
   castep
   cp2k
   crystal
   demon
   demonnano
   dftb
   dmol
   espresso
   exciting
   FHI-aims
   fleur
   gamess_us
   gaussian
   gromacs
   gulp
   harmonic
   socketio/socketio
   jacapo
   kim
   lammps
   mopac
   nwchem
   octopus
   onetep
   openmx
   orca
   plumed
   psi4
   qchem
   siesta
   turbomole
   tersoff
   vasp
   qmmm
   checkpointing
   fd
   loggingcalc
   mixing
   dftd3
   others
   test
   ace
