import sys

from ase.build import molecule
from ase.calculators.siesta import Siesta
from ase.calculators.socketio import SocketIOCalculator
from ase.optimize import BFGS

unixsocket = 'siesta'

fdf_arguments = {
    'MD.TypeOfRun': 'Master',
    'Master.code': 'i-pi',
    'Master.interface': 'socket',
    'Master.address': unixsocket,
    'Master.socketType': 'unix',
}

# To connect through INET socket instead, use:
#   fdf_arguments['Master.port'] = port
#   fdf_arguments['Master.socketType'] = 'inet'
# Optional, for networking:
#   fdf_arguments['Master.address'] = <hostname or IP address>

atoms = molecule('H2O', vacuum=3.0)
atoms.rattle(stdev=0.1)

siesta = Siesta(fdf_arguments=fdf_arguments)
opt = BFGS(atoms, trajectory='opt.siesta.traj', logfile='opt.siesta.log')

with SocketIOCalculator(siesta, log=sys.stdout, unixsocket=unixsocket) as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)

# Note: Siesta does not exit cleanly - expect nonzero exit codes.
