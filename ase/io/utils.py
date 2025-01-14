from itertools import islice
from math import sqrt
from typing import IO

import numpy as np

from ase.data import atomic_numbers, covalent_radii
from ase.data.colors import jmol_colors
from ase.io.formats import string2index
from ase.utils import rotate


class PlottingVariables:
    # removed writer - self
    def __init__(self, atoms, rotation='', show_unit_cell=2,
                 radii=None, bbox=None, colors=None, scale=20,
                 maxwidth=500, extra_offset=(0., 0.)):
        self.numbers = atoms.get_atomic_numbers()
        self.colors = colors
        if colors is None:
            ncolors = len(jmol_colors)
            self.colors = jmol_colors[self.numbers.clip(max=ncolors - 1)]

        if radii is None:
            radii = covalent_radii[self.numbers]
        elif isinstance(radii, float):
            radii = covalent_radii[self.numbers] * radii
        else:
            radii = np.array(radii)

        natoms = len(atoms)

        if isinstance(rotation, str):
            rotation = rotate(rotation)

        cell = atoms.get_cell()
        disp = atoms.get_celldisp().flatten()

        if show_unit_cell > 0:
            L, T, D = cell_to_lines(self, cell)
            cell_vertices = np.empty((2, 2, 2, 3))
            for c1 in range(2):
                for c2 in range(2):
                    for c3 in range(2):
                        cell_vertices[c1, c2, c3] = np.dot([c1, c2, c3],
                                                           cell) + disp
            cell_vertices.shape = (8, 3)
            cell_vertices = np.dot(cell_vertices, rotation)
        else:
            L = np.empty((0, 3))
            T = None
            D = None
            cell_vertices = None

        nlines = len(L)

        positions = np.empty((natoms + nlines, 3))
        R = atoms.get_positions()
        positions[:natoms] = R
        positions[natoms:] = L

        r2 = radii**2
        for n in range(nlines):
            d = D[T[n]]
            if ((((R - L[n] - d)**2).sum(1) < r2) &
                    (((R - L[n] + d)**2).sum(1) < r2)).any():
                T[n] = -1

        positions = np.dot(positions, rotation)
        R = positions[:natoms]

        if bbox is None:
            X1 = (R - radii[:, None]).min(0)
            X2 = (R + radii[:, None]).max(0)
            if show_unit_cell == 2:
                X1 = np.minimum(X1, cell_vertices.min(0))
                X2 = np.maximum(X2, cell_vertices.max(0))
            M = (X1 + X2) / 2
            S = 1.05 * (X2 - X1)
            w = scale * S[0]
            if w > maxwidth:
                w = maxwidth
                scale = w / S[0]
            h = scale * S[1]
            offset = np.array([scale * M[0] - w / 2, scale * M[1] - h / 2, 0])
        else:
            w = (bbox[2] - bbox[0]) * scale
            h = (bbox[3] - bbox[1]) * scale
            offset = np.array([bbox[0], bbox[1], 0]) * scale

        offset[0] = offset[0] - extra_offset[0]
        offset[1] = offset[1] - extra_offset[1]
        self.w = w + extra_offset[0]
        self.h = h + extra_offset[1]

        positions *= scale
        positions -= offset

        if nlines > 0:
            D = np.dot(D, rotation)[:, :2] * scale

        if cell_vertices is not None:
            cell_vertices *= scale
            cell_vertices -= offset

        cell = np.dot(cell, rotation)
        cell *= scale

        self.cell = cell
        self.positions = positions
        self.D = D
        self.T = T
        self.cell_vertices = cell_vertices
        self.natoms = natoms
        self.d = 2 * scale * radii
        self.constraints = atoms.constraints

        # extension for partial occupancies
        self.frac_occ = False
        self.tags = None
        self.occs = None

        try:
            self.occs = atoms.info['occupancy']
            self.tags = atoms.get_tags()
            self.frac_occ = True
        except KeyError:
            pass


def cell_to_lines(writer, cell):
    # XXX this needs to be updated for cell vectors that are zero.
    # Cannot read the code though!  (What are T and D? nn?)
    nlines = 0
    nsegments = []
    for c in range(3):
        d = sqrt((cell[c]**2).sum())
        n = max(2, int(d / 0.3))
        nsegments.append(n)
        nlines += 4 * n

    positions = np.empty((nlines, 3))
    T = np.empty(nlines, int)
    D = np.zeros((3, 3))

    n1 = 0
    for c in range(3):
        n = nsegments[c]
        dd = cell[c] / (4 * n - 2)
        D[c] = dd
        P = np.arange(1, 4 * n + 1, 4)[:, None] * dd
        T[n1:] = c
        for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            n2 = n1 + n
            positions[n1:n2] = P + i * cell[c - 2] + j * cell[c - 1]
            n1 = n2

    return positions, T, D


def make_patch_list(writer):
    from matplotlib.patches import Circle, PathPatch, Wedge
    from matplotlib.path import Path

    indices = writer.positions[:, 2].argsort()
    patch_list = []
    for a in indices:
        xy = writer.positions[a, :2]
        if a < writer.natoms:
            r = writer.d[a] / 2
            if writer.frac_occ:
                site_occ = writer.occs[str(writer.tags[a])]
                # first an empty circle if a site is not fully occupied
                if (np.sum([v for v in site_occ.values()])) < 1.0:
                    # fill with white
                    fill = '#ffffff'
                    patch = Circle(xy, r, facecolor=fill,
                                   edgecolor='black')
                    patch_list.append(patch)

                start = 0
                # start with the dominant species
                for sym, occ in sorted(site_occ.items(),
                                       key=lambda x: x[1],
                                       reverse=True):
                    if np.round(occ, decimals=4) == 1.0:
                        patch = Circle(xy, r, facecolor=writer.colors[a],
                                       edgecolor='black')
                        patch_list.append(patch)
                    else:
                        # jmol colors for the moment
                        extent = 360. * occ
                        patch = Wedge(
                            xy, r, start, start + extent,
                            facecolor=jmol_colors[atomic_numbers[sym]],
                            edgecolor='black')
                        patch_list.append(patch)
                        start += extent

            else:
                if ((xy[1] + r > 0) and (xy[1] - r < writer.h) and
                        (xy[0] + r > 0) and (xy[0] - r < writer.w)):
                    patch = Circle(xy, r, facecolor=writer.colors[a],
                                   edgecolor='black')
                    patch_list.append(patch)
        else:
            a -= writer.natoms
            c = writer.T[a]
            if c != -1:
                hxy = writer.D[c]
                patch = PathPatch(Path((xy + hxy, xy - hxy)))
                patch_list.append(patch)
    return patch_list


class ImageChunk:
    """Base Class for a file chunk which contains enough information to
    reconstruct an atoms object."""

    def build(self, **kwargs):
        """Construct the atoms object from the stored information,
        and return it"""


class ImageIterator:
    """Iterate over chunks, to return the corresponding Atoms objects.
    Will only build the atoms objects which corresponds to the requested
    indices when called.
    Assumes ``ichunks`` is in iterator, which returns ``ImageChunk``
    type objects. See extxyz.py:iread_xyz as an example.
    """

    def __init__(self, ichunks):
        self.ichunks = ichunks

    def __call__(self, fd: IO, index=None, **kwargs):
        if isinstance(index, str):
            index = string2index(index)

        if index is None or index == ':':
            index = slice(None, None, None)

        if not isinstance(index, (slice, str)):
            index = slice(index, (index + 1) or None)

        for chunk in self._getslice(fd, index):
            yield chunk.build(**kwargs)

    def _getslice(self, fd: IO, indices: slice):
        try:
            iterator = islice(self.ichunks(fd),
                              indices.start, indices.stop,
                              indices.step)
        except ValueError:
            # Negative indices.  Go through the whole thing to get the length,
            # which allows us to evaluate the slice, and then read it again
            if not hasattr(fd, 'seekable') or not fd.seekable():
                raise ValueError('Negative indices only supported for '
                                 'seekable streams')

            startpos = fd.tell()
            nchunks = 0
            for _ in self.ichunks(fd):
                nchunks += 1
            fd.seek(startpos)
            indices_tuple = indices.indices(nchunks)
            iterator = islice(self.ichunks(fd), *indices_tuple)
        return iterator


def verify_cell_for_export(cell, check_orthorhombric=True):
    """Function to verify if the cell size is defined and if the cell is

    Parameters:

    cell: cell object
        cell to be checked.

    check_orthorhombric: bool
        If True, check if the cell is orthorhombric, raise an ``ValueError`` if
        the cell is orthorhombric. If False, doesn't check if the cell is
        orthorhombric.

    Raise a ``ValueError`` if the cell if not suitable for export to mustem xtl
    file or prismatic/computem xyz format:
        - if cell is not orthorhombic (only when check_orthorhombric=True)
        - if cell size is not defined
    """

    if check_orthorhombric and not cell.orthorhombic:
        raise ValueError('To export to this format, the cell needs to be '
                         'orthorhombic.')
    if cell.rank < 3:
        raise ValueError('To export to this format, the cell size needs '
                         'to be set: current cell is {}.'.format(cell))


def verify_dictionary(atoms, dictionary, dictionary_name):
    """
    Verify a dictionary have a key for each symbol present in the atoms object.

    Parameters:

    dictionary: dict
        Dictionary to be checked.


    dictionary_name: dict
        Name of the dictionary to be displayed in the error message.

    cell: cell object
        cell to be checked.


    Raise a ``ValueError`` if the key doesn't match the atoms present in the
    cell.
    """
    # Check if we have enough key
    for key in set(atoms.symbols):
        if key not in dictionary:
            raise ValueError('Missing the {} key in the `{}` dictionary.'
                             ''.format(key, dictionary_name))


def segment_list(data, segment_size):
    """Segments a list into sublists of a specified size."""
    return [data[i:i + segment_size] for i in range(0, len(data), segment_size)]
