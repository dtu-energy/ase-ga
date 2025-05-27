# creates: diffusion-I.png, diffusion-T.png, diffusion-F.png
# creates: diffusion-barrier.png
import runpy

from ase.io import read, write
from ase.mep import NEBTools

runpy.run_path('diffusion1.py')
runpy.run_path('diffusion2.py')
runpy.run_path('diffusion4.py')
runpy.run_path('diffusion5.py')

images = read('neb.traj@-5:')
for name, a in zip('ITF', images[::2]):
    cell = a.get_cell()
    del a.constraints
    a = a * (2, 2, 1)
    a.set_cell(cell)
    renderer = write(
        'diffusion-%s.pov' % name,
        a,
        povray_settings=dict(transparent=False, display=False),
    )
    renderer.render()

nebtools = NEBTools(images)
assert abs(nebtools.get_barrier()[0] - 0.374) < 1e-3
