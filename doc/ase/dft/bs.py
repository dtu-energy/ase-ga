# creates: cu.png
from ase.build import bulk
from ase.calculators.test import FreeElectrons

a = bulk('Cu')
a.calc = FreeElectrons(nvalence=1, kpts={'path': 'GXWLGK', 'npoints': 200})
a.get_potential_energy()
bs = a.calc.band_structure()
bs.plot(emin=0, emax=20, filename='cu.png')
