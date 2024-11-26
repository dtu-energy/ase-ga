import matplotlib.pyplot as plt
import numpy as np
from gpaw import GPAW

from ase import Atom, Atoms

ase = """\
 H   HH HHH
H H H   H
HHH  H  HH
H H   H H
H H HH  HHH"""

d = 1.2

logo = Atoms()
for i, line in enumerate(ase.split('\n')):
    for j, c in enumerate(line):
        if c == 'H':
            logo.append(Atom('H', [d * j, d * i, 0]))
logo.set_cell((15, 15, 2))
logo.center()

calc = GPAW()
logo.calc = calc
e = logo.get_potential_energy()
calc.write('logo2.gpw')

print(calc.density.nt_sg.shape)
n = calc.density.nt_sg[0, :, :, 10]
# 1c4e63
c0 = np.array([19, 63, 82.0]).reshape((3, 1, 1)) / 255
c1 = np.array([1.0, 1, 0]).reshape((3, 1, 1))
a = c0 + n / n.max() * (c1 - c0)
print(a.shape)

i = plt.imshow(a.T, aspect=True)
i.write_png('ase.png')
plt.axis('off')
plt.savefig('ase2.png', dpi=200)
