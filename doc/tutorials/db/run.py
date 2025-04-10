# creates: cu1o.png, cu2o.png, cu3o.png
import os
import runpy

from ase.cli.main import main
from ase.db import connect
from ase.io import read, write

for name in ['bulk.db', 'ads.db', 'refs.db']:
    if os.path.isfile(name):
        os.remove(name)

# Run the tutorial:
for filename in ['bulk.py', 'ads.py', 'refs.py']:
    runpy.run_path(filename)

for cmd in [
    'ase db ads.db ads=clean --insert-into refs.db',
    'ase db ads.db ads=clean --delete --yes',
    'ase db ads.db pbc=FFF --insert-into refs.db',
    'ase db ads.db pbc=FFF --delete --yes',
]:
    main(args=cmd.split()[1:])

runpy.run_path('ea.py')

# Create the figures:
for n in [1, 2, 3]:
    a = read(f'ads.db@Cu{n}O')[0]
    a *= (2, 2, 1)
    renderer = write(f'cu{n}o.pov', a, rotation='-80x')
    renderer.render()

# A bit of testing:
row = connect('ads.db').get(surf='Pt', layers=3, ads='O')
assert abs(row.ea - -4.724) < 0.001
assert abs(row.height - 1.706) < 0.001
