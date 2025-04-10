# creates: surface.png
import runpy

from ase.io import read, write

runpy.run_path('N2Cu.py')
image = read('N2Cu.traj@-1')
write('surface.pov', image, povray_settings=dict(transparent=False)).render()
