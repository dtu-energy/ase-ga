import matplotlib.pyplot as plt
from gpaw import GPAW

from ase.dft.dos import DOS

calc = GPAW('bulk.Ag.gpw')

dos = DOS(calc, npts=800, width=0)
energies = dos.get_energies()
weights = dos.get_dos()

ax = plt.gca()
ax.plot(energies, weights)
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('DOS [1/eV]')
plt.savefig('dos.png')
plt.show()
