# creates: saxs.png, xrd.png

import numpy as np

from ase.cluster.cubic import FaceCenteredCubic
from ase.utils.xrdebye import XrDebye

# create nanoparticle with approx. 2 nm diameter
atoms = FaceCenteredCubic(
    'Ag', [(1, 0, 0), (1, 1, 0), (1, 1, 1)], [6, 8, 8], 4.09
)
# setup for desired wavelength
xrd = XrDebye(atoms=atoms, wavelength=0.50523)
# calculate and plot diffraction pattern
xrd.calc_pattern(x=np.arange(15, 30, 0.1), mode='XRD')
xrd.plot_pattern('xrd.png')
# calculate and plot samll-angle scattering
xrd.calc_pattern(x=np.logspace(-2, -0.3, 50), mode='SAXS')
xrd.plot_pattern('saxs.png')
