# creates: matplotlib_plot_atoms1.png, matplotlib_plot_atoms2.png
# creates: matplotlib_plot_atoms3.png

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from ase.lattice.cubic import FaceCenteredCubic
from ase.spacegroup import crystal
from ase.visualize.plot import plot_atoms

slab = FaceCenteredCubic('Au', size=(2, 2, 2))
fig, ax = plt.subplots()
plot_atoms(slab, ax, radii=0.3, rotation=('90x,45y,0z'))
fig.savefig('matplotlib_plot_atoms1.png')

slab = FaceCenteredCubic('Au', size=(2, 2, 2))
fig, axarr = plt.subplots(1, 4, figsize=(15, 5))
plot_atoms(slab, axarr[0], radii=0.3, rotation=('0x,0y,0z'))
plot_atoms(
    slab, axarr[1], scale=0.7, offset=(3, 4), radii=0.3, rotation=('0x,0y,0z')
)
plot_atoms(slab, axarr[2], radii=0.3, rotation=('45x,45y,0z'))
plot_atoms(slab, axarr[3], radii=0.3, rotation=('0x,0y,0z'))
axarr[0].set_title('No rotation')
axarr[1].set_xlabel(r'X-axis, [$\mathrm{\AA}$]')
axarr[1].set_ylabel(r'Y-axis, [$\mathrm{\AA}$]')
axarr[2].set_axis_off()
axarr[3].set_xlim(2, 6)
axarr[3].set_ylim(2, 6)
fig.savefig('matplotlib_plot_atoms2.png')

stem_image = mpimg.imread('stem_image.jpg')
atom_pos = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5), (0.5, 0.5, 0.0)]
srtio3 = crystal(
    ['Sr', 'Ti', 'O'], atom_pos, spacegroup=221, cellpar=3.905, size=(3, 3, 3)
)
fig, ax = plt.subplots()
ax.imshow(stem_image, cmap='gray')
plot_atoms(
    srtio3, ax, radii=0.3, scale=6.3, offset=(47, 54), rotation=('90x,45y,56z')
)
ax.set_xlim(0, stem_image.shape[0])
ax.set_ylim(0, stem_image.shape[1])
ax.set_axis_off()
fig.tight_layout()
fig.savefig('matplotlib_plot_atoms3.png')
