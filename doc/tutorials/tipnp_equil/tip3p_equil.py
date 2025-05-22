import numpy as np

import ase.units as units
from ase import Atoms
from ase.calculators.tip3p import TIP3P, angleHOH, rOH
from ase.constraints import FixBondLengths
from ase.io.trajectory import Trajectory
from ase.md import Langevin

# Set up water box at 20 deg C density
x = angleHOH * np.pi / 180 / 2
pos = [
    [0, 0, 0],
    [0, rOH * np.cos(x), rOH * np.sin(x)],
    [0, rOH * np.cos(x), -rOH * np.sin(x)],
]
atoms = Atoms('OH2', positions=pos)

vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24)) ** (1 / 3.0)
atoms.set_cell((vol, vol, vol))
atoms.center()

atoms = atoms.repeat((3, 3, 3))
atoms.set_pbc(True)

# RATTLE-type constraints on O-H1, O-H2, H1-H2.
atoms.constraints = FixBondLengths(
    [(3 * i + j, 3 * i + (j + 1) % 3) for i in range(3**3) for j in [0, 1, 2]]
)

tag = 'tip3p_27mol_equil'
atoms.calc = TIP3P(rc=4.5)
md = Langevin(
    atoms,
    1 * units.fs,
    temperature=300 * units.kB,
    friction=0.01,
    logfile=tag + '.log',
)

traj = Trajectory(tag + '.traj', 'w', atoms)
md.attach(traj.write, interval=1)
md.run(4000)

# Repeat box and equilibrate further.
tag = 'tip3p_216mol_equil'
atoms.set_constraint()  # repeat not compatible with FixBondLengths currently.
atoms = atoms.repeat((2, 2, 2))
atoms.constraints = FixBondLengths(
    [
        (3 * i + j, 3 * i + (j + 1) % 3)
        for i in range(len(atoms) / 3)
        for j in [0, 1, 2]
    ]
)
atoms.calc = TIP3P(rc=7.0)
md = Langevin(
    atoms,
    2 * units.fs,
    temperature=300 * units.kB,
    friction=0.01,
    logfile=tag + '.log',
)

traj = Trajectory(tag + '.traj', 'w', atoms)
md.attach(traj.write, interval=1)
md.run(2000)
