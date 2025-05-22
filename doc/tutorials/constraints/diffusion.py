# creates:  diffusion-path.png
import runpy

from ase.io import read, write

runpy.run_path('diffusion4.py')
images = [read('mep%d.traj' % i) for i in range(5)]
a = images[0] + images[1] + images[2] + images[3] + images[4]
del a.constraints
a *= (2, 1, 1)
a.set_cell(images[0].get_cell())
renderer = write(
    'diffusion-path.pov',
    a,
    rotation='-90x',
    povray_settings=dict(transparent=False),
)
renderer.render()
