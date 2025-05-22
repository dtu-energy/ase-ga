from gpaw import GPAW, FermiDirac
from gpaw.analyse.overlap import Overlap
from gpaw.cluster import Cluster
from gpaw.lrtddft import LrTDDFT

from ase.vibrations.resonant_raman import ResonantRamanCalculator

h = 0.25
atoms = Cluster('relaxed.traj')
atoms.minimal_box(3.5, h=h)

# relax the molecule
calc = GPAW(
    h=h,
    occupations=FermiDirac(width=0.1),
    symmetry={'point_group': False},
    nbands=10,
    convergence={'eigenstates': 1.0e-5, 'bands': 4},
)
atoms.calc = calc

# use only the 4 converged states for linear response calculation
rmc = ResonantRamanCalculator(
    atoms,
    LrTDDFT,
    name='rroverlap',
    exkwargs={'restrict': {'jend': 3}},
    overlap=lambda x, y: Overlap(x).pseudo(y)[0],
)
rmc.run()
