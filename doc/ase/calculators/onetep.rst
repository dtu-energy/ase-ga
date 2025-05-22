.. module:: ase.calculators.onetep

======
ONETEP
======

Introduction
============

ONETEP_ is a linear-scaling density functional theory code which exploit the
near-sightness of the electronic density. It uses a set of atom-centered local 
orbitals (denoted NGWFs) which are optimised in situ to enable calculations 
with a minimal number of orbitals.

This interface makes it possible to use ONETEP as a calculator in ASE.
You need to have a copy of the ONETEP code (and an appropriate license) to use
this interface.

.. _ONETEP: http://www.onetep.org

Environment variables
=====================

The environment variable :envvar:`ASE_ONETEP_COMMAND` must hold the command
to invoke the ONETEP calculation. The variable must be a string with a link
to the ONETEP binary, and any other specific settings required for your
environment (srun, mpirun, ...)

You can setup this environment variable in your shell configuration file:

.. highlight:: bash

::

  $ export ASE_ONETEP_COMMAND="export OMP_NUM_THREADS=4; mpirun -n 6 ~/onetep/bin/onetep.arch"


.. highlight:: python

Or within python itself:

  >>> environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=4; mpirun -n 6 ~/onetep/bin/onetep.arch"

ASE will automatically redirect stdout and stderr to the appropriate
files, namely "$LABEL.out" and "$LABEL.err" where label is the name
used for your ONETEP calculations

Pseudopotentials
================

ONETEP accepts PAW datasets and NCP pseudpotentials with formats 
USP and recpot. Pseudopotentials are passed directly to the Onetep calculator
as a dictionnary definition. If no pseudopotentials are passed ASE will
try to guess the files based on the element used and the pseudo_path variable. ::

    # Explicitly providing each path
    calc = Onetep(pseudopotentials = {'H': '/path/to/pseudos/H.usp', 'O': '/path/to/pseudos/O.usp'})
    # Using pseudo_path
    calc = Onetep(pseudo_path = '/path/to/pseudos', pseudopotentials = {'H': 'H.usp', 'O': 'O.usp'})
    # ASE will try to guess them
    calc = Onetep(pseudo_path = '/path/to/pseudos')

For ASE to correctly guess the pseudopotentials use a pseudo_path that contains only one kind
of pseudopotential per element.

.. highlight:: python

ONETEP Calculator
=================

.. autoclass:: ase.calculators.onetep.Onetep

Simple calculations can be setup calling the Onetep calculator without any parameters,
in this case ONETEP default parameters will be used. For more complex cases using the
keywords parameters is necessary. The 'keywords' parameters is a dictionnary in which
each key is a string that represent a ONETEP keywords.

Examples
========

Here is an example of setting up a calculation on a water molecule: ::

    from ase.build import molecule
    from ase.calculators.onetep import Onetep
    from os import environ

    # water molecule from ASE database, centered in a ~ 24 Å box
    wat = molecule('H2O')
    wat.center(12)
    environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=8; mpirun -n 2 ~/onetep/bin/onetep.arch"
    # Ouput will be in "water.out"
    calc = Onetep(label = 'water', xc = 'PBE', paw = True, pseudo_path = '/path/to/pseudos')
    wat.calc = calc
    wat.get_potential_energy()

.. highlight:: python

Another more complex example on Pt13: ::

    from os import environ

    import numpy as np

    from ase.build import molecule
    from ase.calculators.onetep import Onetep
    from ase.cluster import Octahedron
    from ase.optimize.sciopt import SciPyFminBFGS
    # Pt13 from ase.cluster
    nano = Octahedron('Pt', 3, 1)
    nano.set_cell(np.eye(3)*12)
    nano.center()

    label = 'pt13'

    environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=8; mpirun -n 8 ~/onetep/bin/onetep.arch"

    # ONETEP default are atomic units, one can specify 'cutoff_energy' : '600 eV' if needed.
    keywords = {
        'xc' : 'rpbe',
        'do_properties' : True,
        'cutoff_energy' : 35,
        'output_detail': 'verbose',
        'elec_energy_tol': 1.0e-5/len(atoms),
    }

    # Ouput will be in "pt13.out", 
    # append = True will not overwrite file at each step
    calc = Onetep(
        label = label,
        edft = True,
        append = True,
        pseudo_path = '/path/to/pseudos', 
        keywords = keywords)

    nanoparticle.calc = calc

    opt = SciPyFminBFGS(atoms = nano, trajectory = label + ".traj", logfile = label + ".log")
    opt.run(fmax=0.01)

.. highlight:: python

Here is an example of setting up an EELS and LDOS calculations on an N-substituted graphene sheet,
demonstrating several more advanced functionalities (eg tags, species groups, and overrides to
pseudopotentials and atomic solver strings): ::

    import numpy as np

    from ase.build import graphene_nanoribbon
    from ase.calculators.onetep import Onetep
    from ase.io import write
    from numpy.linalg import norm
    from numpy.random import choice

    sheet = graphene_nanoribbon(10, 10, type='zigzag', vacuum = 10)

    # Get all distances to center of mass
    com = sheet.get_center_of_mass()
    distances_to_com = norm(sheet.positions - com, axis = 1)

    # Find atoms close to com and change one randomly to N
    p, = np.where(distances_to_com < 5)
    to_nitro = choice(p)
    sheet[to_nitro].symbol = 'N'

    shell_rad = np.array([1.5, 2.5, 3.0, 4.0, 4.5])

    tags = np.zeros(len(sheet), dtype=np.int32)

    # We want to tag atoms that are close to the introduced nitrogen
    for idx, rad in enumerate(reversed(shell_rad)):
        # All distances N-C
        dist = norm(sheet[to_nitro].position - sheet.get_positions(), axis = 1)
        # Which ones are closest to rad?
        p, = np.where(dist < rad)
        # Cannot be the nitrogen itself
        p = p[p != to_nitro]
        # Tags them
        tags[p] = len(shell_rad) - idx

    sheet.set_tags(tags)

    tags = ['' if i == 0 else i for i in tags]

    species = np.unique(np.char.add(sheet.get_chemical_symbols(), tags))

    keywords = {
        'species_core_wf' : ['N /path/to/pseudo/corehole.abinit'],
        'species_solver' : ['N SOLVE conf=1s1 2p4'],
        'pseudo_path': '/Users/tomdm/PseudoPotentials/SSSP_1.2.1',
        'xc' : 'PBE',
        'paw': True,
        'do_properties': True,
        'cutoff_energy' : '500 eV',
        'species_ldos_groups': species,
        'task' : 'GeometryOptimization'
    }

    calc = Onetep(
        label = 'N_doped_graphene_001',
        keywords = keywords
    )

    # Checking the input before running the calculation
    write('to_check.dat', sheet, format='onetep-in', keywords = keywords)

    sheet.calc = calc
    # Will actually run the geometry optimisation
    # using ONETEP internal BFGS
    sheet.get_potential_energy()

.. highlight:: python

Quickly restart with solvation effect using the soft sphere model ::

    from ase.io import read
    from ase.io.onetep import get_onetep_keywords

    # Read from the previous run...
    optimized_sheet = read("N_doped_graphene_001.out")

    # Function to retrieve keywords dict from input file...
    keywords = get_onetep_keywords('N_doped_graphene_001.dat')

    # We add solvation keywords
    keywords.update(
        {
        'is_implicit_solvent': True,
        'is_include_apolar': True,
        'is_smeared_ion_rep': True,
        'is_dielectric_model': 'fix_cavity',
        'is_dielectric_function' : 'soft_sphere'
        }
    )

    optimized_sheet.calc = Onetep(...)

    ...


