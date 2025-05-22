# creates: spacegroup-al.png spacegroup-fe.png spacegroup-rutile.png spacegroup-cosb3.png spacegroup-mg.png spacegroup-skutterudite.png spacegroup-diamond.png spacegroup-nacl.png
# flake8: noqa
import runpy

import ase.io

for name in ['al', 'mg', 'fe', 'diamond', 'nacl', 'rutile', 'skutterudite']:
    py = f'spacegroup-{name}.py'
    dct = runpy.run_path(py)
    atoms = dct[name]
    renderer = ase.io.write(
        'spacegroup-%s.pov' % name,
        atoms,
        rotation='10x,-10y',
        povray_settings=dict(
            transparent=False,
            # canvas_width=128,
            # celllinewidth=0.02,
            celllinewidth=0.05,
        ),
    )
    renderer.render()

runpy.run_path('spacegroup-cosb3.py')
