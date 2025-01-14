# creates: zno.png
import numpy as np

from ase.phasediagram import solvated
from ase.pourbaix import Pourbaix

refs = solvated('Zn')
refs += [('Zn', 0.0), ('ZnO', -3.323), ('ZnO2(aq)', -2.921)]
pbx = Pourbaix('ZnO', dict(refs))
energy, phase = pbx.get_pourbaix_energy(1.0, 9.0, verbose=True)
print(type(phase))
print(phase.get_free_energy(1.0, 9.0))
pbx.get_pourbaix_energy(0.0, 10.0, verbose=True)

diagram = pbx.diagram(U=np.linspace(-2, 2, 100), pH=np.linspace(0, 14, 100))

diagram.plot(
    show=False,
    include_text=True,
    filename='zno.png'
)
