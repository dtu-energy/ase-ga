# creates: a1.png, a2.png, a3.png, cnt1.png, cnt2.png, gnr1.png, gnr2.png
from ase.build import bulk, graphene_nanoribbon, nanotube
from ase.io import write

for i, a in enumerate(
    [bulk('Cu', 'fcc', a=3.6),
     bulk('Cu', 'fcc', a=3.6, orthorhombic=True),
     bulk('Cu', 'fcc', a=3.6, cubic=True)]):
    write('a%d.pov' % (i + 1), a).render()

cnt1 = nanotube(6, 0, length=4, vacuum=2.5)
cnt1.rotate('x', 'z', rotate_cell=True)
cnt2 = nanotube(3, 3, length=6, bond=1.4, symbol='Si', vacuum=2.5)
cnt2.rotate('x', 'z', rotate_cell=True)

for i, a in enumerate([cnt1, cnt2]):
    write('cnt%d.pov' % (i + 1), a).render()

gnr1 = graphene_nanoribbon(3, 4, type='armchair', saturated=True, vacuum=2.5)
gnr2 = graphene_nanoribbon(2, 6, type='zigzag', saturated=True,
                           C_H=1.1, C_C=1.4, vacuum=3.0,
                           magnetic=True, initial_mag=1.12)

for i, a in enumerate([gnr1, gnr2]):
    write('gnr%d.pov' % (i + 1), a,
          rotation='90x').render()
