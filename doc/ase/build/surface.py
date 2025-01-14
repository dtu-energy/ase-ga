# creates: fcc100.png, fcc110.png, bcc100.png, fcc111.png, bcc110.png
# creates: bcc111.png, hcp0001.png, fcc111o.png, fcc211o.png, bcc110o.png
# creates: bcc111o.png, hcp0001o.png, ontop-site.png, hollow-site.png
# creates: fcc-site.png, hcp-site.png, bridge-site.png, diamond100.png
# creates: diamond111.png, hcp10m10.png, mx2.png, fcc111_root.png
# creates: graphene.png

import ase.build as surface
from ase import Atoms
from ase.build import fcc111, root_surface
from ase.io import write

surfaces = ['fcc100', 'fcc110', 'bcc100', 'hcp10m10', 'diamond100',
            'fcc111', 'bcc110', 'bcc111', 'hcp0001', 'diamond111', 'fcc211',
            'mx2', 'graphene']

symbols = {
    'fcc': 'Cu',
    'bcc': 'Fe',
    'hcp': 'Ru',
    'dia': 'C',
    'mx2': 'MoS2',
    'gra': 'C2'}
radii = {
    'fcc': 1.1,
    'bcc': 1.06,
    'hcp': 1.08,
    'dia': 0.5,
    'mx2': 1.0,
    'gra': 1.0}
adsorbates = {'ontop': 'H', 'hollow': 'O', 'fcc': 'N', 'hcp': 'C',
              'bridge': 'F'}


def save(name, slab):
    print('save %s' % name)
    write(name + '.png', slab, radii=radii[name[:3]],
          scale=10)


for name in surfaces:
    f = getattr(surface, name)
    for kwargs in [{}, {'orthogonal': False}, {'orthogonal': True}]:
        print(name, kwargs)
        try:
            slab = f(symbols[name[:3]], size=(3, 4, 5), vacuum=4, **kwargs)
        except (TypeError, NotImplementedError):
            continue
        try:
            for site in slab.info['adsorbate_info']['sites']:
                if site.endswith('bridge'):
                    h = 1.5
                else:
                    h = 1.2
                surface.add_adsorbate(slab, adsorbates.get(site, 'F'), h, site)
        except KeyError:
            pass
        if kwargs.get('orthogonal', None):
            name += 'o'
        save(name, slab)

for site, symbol in adsorbates.items():
    write('%s-site.png' % site, Atoms(symbol), radii=1.08, scale=10)

fcc111_primitive = fcc111('Ag', (1, 1, 3))
fcc111_root = root_surface(fcc111_primitive, 27)
save('fcc111_root', fcc111_root)
