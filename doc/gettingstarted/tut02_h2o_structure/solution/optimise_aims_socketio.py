import os

from ase import Atoms
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from ase.optimize import BFGS

os.environ['ASE_AIMS_COMMAND'] = 'aims.x'
os.environ['AIMS_SPECIES_DIR'] = '/home/myname/FHIaims/species_defaults/light'

atoms = Atoms('HOH', positions=[[0, 0, -1], [0, 1, 0], [0, 0, 1]])
opt = BFGS(atoms, trajectory='opt-aims-socketio.traj')

aims = Aims(
    xc='LDA', compute_forces=True, use_pimd_wrapper=('UNIX:mysocket', 31415)
)

with SocketIOCalculator(aims, unixsocket='mysocket') as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)
