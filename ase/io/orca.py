import re
from io import StringIO
from pathlib import Path
from typing import List, Optional

import numpy as np

from ase.io import read
from ase.units import Bohr, Hartree
from ase.utils import reader, writer

# Made from NWChem interface


@reader
def read_geom_orcainp(fd):
    """Method to read geometry from an ORCA input file."""
    lines = fd.readlines()

    # Find geometry region of input file.
    stopline = 0
    for index, line in enumerate(lines):
        if line[1:].startswith('xyz '):
            startline = index + 1
            stopline = -1
        elif (line.startswith('end') and stopline == -1):
            stopline = index
        elif (line.startswith('*') and stopline == -1):
            stopline = index
    # Format and send to read_xyz.
    xyz_text = '%i\n' % (stopline - startline)
    xyz_text += ' geometry\n'
    for line in lines[startline:stopline]:
        xyz_text += line
    atoms = read(StringIO(xyz_text), format='xyz')
    atoms.set_cell((0., 0., 0.))  # no unit cell defined

    return atoms


@writer
def write_orca(fd, atoms, params):
    # conventional filename: '<name>.inp'
    fd.write(f"! {params['orcasimpleinput']} \n")
    fd.write(f"{params['orcablocks']} \n")

    if 'coords' not in params['orcablocks']:
        fd.write('*xyz')
        fd.write(" %d" % params['charge'])
        fd.write(" %d \n" % params['mult'])
        for atom in atoms:
            if atom.tag == 71:  # 71 is ascii G (Ghost)
                symbol = atom.symbol + ' : '
            else:
                symbol = atom.symbol + '   '
            fd.write(
                symbol
                + str(atom.position[0])
                + " "
                + str(atom.position[1])
                + " "
                + str(atom.position[2])
                + "\n"
            )
        fd.write('*\n')


def read_charge(lines: List[str]) -> Optional[float]:
    """Read sum of atomic charges."""
    charge = None
    for line in lines:
        if 'Sum of atomic charges' in line:
            charge = float(line.split()[-1])
    return charge


def read_energy(lines: List[str]) -> Optional[float]:
    """Read energy."""
    energy = None
    for line in lines:
        if 'FINAL SINGLE POINT ENERGY' in line:
            if "Wavefunction not fully converged" in line:
                energy = float('nan')
            else:
                energy = float(line.split()[-1])
    if energy is not None:
        return energy * Hartree
    return energy


def read_center_of_mass(lines: List[str]) -> Optional[np.ndarray]:
    """ Scan through text for the center of mass """
    # Example:
    # 'The origin for moment calculation is the CENTER OF MASS  =
    # ( 0.002150, -0.296255  0.086315)'
    # Note the missing comma in the output
    com = None
    for line in lines:
        if 'The origin for moment calculation is the CENTER OF MASS' in line:
            line = re.sub(r'[(),]', '', line)
            com = np.array([float(_) for _ in line.split()[-3:]])
    if com is not None:
        return com * Bohr  # return the last match
    return com


def read_dipole(lines: List[str]) -> Optional[np.ndarray]:
    """Read dipole moment.

    Note that the read dipole moment is for the COM frame of reference.
    """
    dipole = None
    for line in lines:
        if 'Total Dipole Moment' in line:
            dipole = np.array([float(_) for _ in line.split()[-3:]])
    if dipole is not None:
        return dipole * Bohr  # Return the last match
    return dipole


@reader
def read_orca_output(fd):
    """ From the ORCA output file: Read Energy and dipole moment
    in the frame of reference of the center of mass "
    """
    lines = fd.readlines()

    energy = read_energy(lines)
    charge = read_charge(lines)
    com = read_center_of_mass(lines)
    dipole = read_dipole(lines)

    results = {}
    results['energy'] = energy
    results['free_energy'] = energy

    if com is not None and dipole is not None:
        dipole = dipole + com * charge
        results['dipole'] = dipole

    return results


@reader
def read_orca_engrad(fd):
    """Read Forces from ORCA .engrad file."""
    getgrad = False
    gradients = []
    tempgrad = []
    for _, line in enumerate(fd):
        if line.find('# The current gradient') >= 0:
            getgrad = True
            gradients = []
            tempgrad = []
            continue
        if getgrad and "#" not in line:
            grad = line.split()[-1]
            tempgrad.append(float(grad))
            if len(tempgrad) == 3:
                gradients.append(tempgrad)
                tempgrad = []
        if '# The at' in line:
            getgrad = False

    forces = -np.array(gradients) * Hartree / Bohr
    return forces


def read_orca_outputs(directory, stdout_path):
    stdout_path = Path(stdout_path)
    results = {}
    results.update(read_orca_output(stdout_path))

    # Does engrad always exist? - No!
    # Will there be other files -No -> We should just take engrad
    # as a direct argument.  Or maybe this function does not even need to
    # exist.
    engrad_path = stdout_path.with_suffix('.engrad')
    if engrad_path.is_file():
        results['forces'] = read_orca_engrad(engrad_path)
    return results
