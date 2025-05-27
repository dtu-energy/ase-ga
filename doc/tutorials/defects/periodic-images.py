# creates: periodic-images-1.svg periodic-images-2.svg

from itertools import product

import matplotlib.pyplot as plt
import numpy as np


class CellFigure:
    def __init__(self, dim):
        """
        Set up a figure for visualizing a cell metric.
        """
        self.fig = plt.figure(figsize=(9, 5))
        self.ax = self.fig.gca()
        self.ax.set_axis_off()
        self.ax.autoscale_view(tight=True)
        self.ax.set_xlim(-2 * dim, 3 * dim)
        self.ax.set_ylim(-dim, 1.5 * dim)
        self.ax.set_aspect('equal')

    def add_cell(
        self,
        cell,
        offset=[0, 0],
        fill_color=None,
        atom=None,
        radius=0.1,
        atom_color='orange',
    ):
        """
        Draw a cell, optionally filled and including an atom.
        """
        xvecs = np.array([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)], dtype=float)
        vectors = []
        for xv in xvecs:
            vectors.append(np.dot(cell, xv + offset))
        vectors = np.array(vectors)
        from matplotlib.patches import Circle, Polygon

        if fill_color is not None:
            self.ax.add_patch(
                Polygon(vectors, closed=True, color=fill_color, alpha=0.4)
            )
        for points in vectors:
            self.ax.plot(vectors.T[0], vectors.T[1], c='k', ls='-')
        if atom:
            pos = np.dot(cell, np.array(atom) + offset)
            self.ax.add_patch(Circle(pos, radius, color=atom_color))
            return pos

    def add_vector(self, cell, xvec, xdir, fill_color='red'):
        """
        Draw an arrow typically symbolizing a neighbor connection.
        """
        from matplotlib.patches import Arrow

        pos = np.dot(cell, xvec)
        dir = np.dot(cell, xdir)
        self.ax.add_patch(Arrow(pos[0], pos[1], dir[0], dir[1], width=0.2))

    def annotate_figure(self, text):
        """
        Add some annotation to the lower left corner of the plot.
        """
        self.ax.text(-2 * dim, 1.35 * dim, text, ha='left', va='center')


# extent of plotted area
dim = 2
# rescaling factor for arrows
rescale = 1.0

# Figure 1
myfig = CellFigure(dim)
prim = np.eye(2)
atompos = [0.5, 0.5]
for i, j in [[-1, 0], [1, 0], [0, -1], [0, 1]]:
    myfig.add_vector(prim, atompos, rescale * np.array([i, j]))
for i, j in product(range(-dim, dim + 1), repeat=2):
    fill_color = 'blue' if i == 0 and j == 0 else None
    myfig.add_cell(prim, [i, j], atom=atompos, fill_color=fill_color)
myfig.annotate_figure('square lattice\n$r_1 = a$, $Z_1=4$')
plt.savefig('periodic-images-1.svg', bbox_inches='tight')

# Figure 2
myfig = CellFigure(dim)
prim = np.array([[2, 0], [0, 0.5]])
for i, j in [[0, -1], [0, 1]]:
    myfig.add_vector(prim, atompos, rescale * np.array([i, j]))
for i, j in product(range(-dim, dim + 1), repeat=2):
    fill_color = 'blue' if i == 0 and j == 0 else None
    myfig.add_cell(prim, [i, j], atom=atompos, fill_color=fill_color)
myfig.annotate_figure(
    'rectangular lattice with a 2:1 aspect ratio\n$r_1 = a/2$, $Z_1=2$'
)
plt.savefig('periodic-images-2.svg', bbox_inches='tight')

# Figure 3
myfig = CellFigure(dim)
prim = np.array([[1, 0.5], [0, np.sqrt(3) / 2]])
prim /= np.linalg.det(prim) ** (1.0 / 2)
s = np.sqrt(2.0 / np.sqrt(3))
for i, j in [[-1, 0], [1, 0], [0, -1], [0, 1], [-1, 1], [1, -1]]:
    myfig.add_vector(prim, atompos, rescale * np.array([i, j]))
for i, j in product(range(-dim, dim + 1), repeat=2):
    fill_color = 'blue' if i == 0 and j == 0 else None
    myfig.add_cell(prim, [i, j], atom=atompos, fill_color=fill_color)
myfig.annotate_figure('hexagonal lattice\n$r_1 = %.3f a$, $Z_1=6$' % s)
plt.savefig('periodic-images-3.svg', bbox_inches='tight')
