# creates: supercell-1.svg supercell-2.svg supercell-3.svg

from math import cos, pi, sin

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def vertices(cell):
    """
    Set up vertices for a cell metric.
    """
    from scipy.spatial import Voronoi

    I = np.indices((3, 3, 3)).reshape((3, 27)) - 1
    G = np.dot(cell, I).T
    vor = Voronoi(G)
    vert1 = []
    for vertices, points in zip(vor.ridge_vertices, vor.ridge_points):
        if -1 not in vertices and 13 in points:
            normal = G[points].sum(0)
            normal /= (normal**2).sum() ** 0.5
            vert1.append((vor.vertices[vertices], normal))
    return vert1


class CellFigure:
    def __init__(self, dim, azim, elev):
        """
        Set up a figure for visualizing a cell metric.
        """
        Axes3D  # silence pyflakes
        self.fig = plt.figure(figsize=(5, 5))
        self.ax = self.fig.add_subplot(projection='3d')
        x = sin(azim)
        y = cos(azim)
        self.view = [x * cos(elev), y * cos(elev), sin(elev)]
        self.ax.set_axis_off()
        self.ax.autoscale_view(tight=True)
        self.ax.set_xlim(0, dim)
        self.ax.set_ylim(0, dim)
        self.ax.set_zlim(0, dim)
        self.ax.view_init(azim=azim / pi * 180, elev=elev / pi * 180)

    def add_cell(self, cell):
        """
        Draw a cell.
        (i.e. edges but no faces)
        """
        vert1 = vertices(cell)
        shift = -vert1[0][0][0]
        for points, normal in vert1:
            if np.dot(normal, self.view) < 0:
                ls = ':'
            else:
                ls = '-'
            x, y, z = np.concatenate([points + shift, points[:1] + shift]).T
            self.ax.plot(x, y, z, c='k', ls=ls)

    def add_primitive_cell(self, cell):
        """
        Draw a primitive unit cell.
        (i.e. a cell with colored faces)
        """
        self.add_cell(cell)
        uc = prim[0][0]
        # plot side faces of unit cell
        X, Y = np.meshgrid([0, uc], [0, uc])
        Z = np.zeros((2, 2)) + uc
        self.ax.plot_surface(
            X, Y, Z, color='blue', alpha=0.5, linewidth=0, zorder=1
        )
        X, Z = np.meshgrid([0, uc], [0, uc])
        Y = np.zeros((2, 2)) + uc
        self.ax.plot_surface(
            X, Y, Z, color='blue', alpha=0.5, linewidth=0, zorder=1
        )
        Y, Z = np.meshgrid([0, uc], [0, uc])
        X = np.zeros((2, 2)) + uc
        self.ax.plot_surface(
            X, Y, Z, color='blue', alpha=0.5, linewidth=0, zorder=1
        )

    def add_atom(self, x0, y0, z0, radius=0.06):
        """
        Draw an atom.
        """
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = x0 + radius * np.outer(np.cos(u), np.sin(v))
        y = y0 + radius * np.outer(np.sin(u), np.sin(v))
        z = z0 + radius * np.outer(np.ones(np.size(u)), np.cos(v))
        self.ax.plot_surface(
            x,
            y,
            z,
            rstride=4,
            cstride=4,
            color='orange',
            linewidth=0.1,
            alpha=0.5,
        )

    def annotate_figure(self, text):
        """
        Add some annotation to the lower left corner of the plot.
        """
        self.ax.text(1.1, 0, -0.2, text, ha='left', va='center')


# extent of plotted area
dim = 0.82
# view angle
azim = 0.75 * pi / 5
elev = 0.5 * pi / 6

# define unit cell and supercell
prim = 1.0 / 3 * np.eye(3)
supr = np.eye(3)

# Figure 1
myfig = CellFigure(dim, azim, elev)
myfig.add_primitive_cell(prim)
myfig.annotate_figure('primitive unit cell')
plt.savefig('supercell-1.svg', bbox_inches='tight')

# Figure 2
myfig = CellFigure(dim, azim, elev)
myfig.add_primitive_cell(prim)
myfig.add_cell(supr)
myfig.annotate_figure('ideal supercell')
plt.savefig('supercell-2.svg', bbox_inches='tight')

# Figure 3
myfig = CellFigure(dim, azim, elev)
myfig.add_cell(supr)
d = 0.08
myfig.add_atom(0.5, 0.5 - d, 0.5)
myfig.add_atom(0.5, 0.5 + d, 0.5)
myfig.annotate_figure('defect supercell')
plt.savefig('supercell-3.svg', bbox_inches='tight')
