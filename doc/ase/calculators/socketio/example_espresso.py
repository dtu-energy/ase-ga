import sys

from ase.build import molecule
from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator
from ase.optimize import BFGS

atoms = molecule('H2O', vacuum=3.0)
atoms.rattle(stdev=0.1)

# Environment-dependent parameters (please configure before running):
pseudopotentials = {'H': 'H.pbe-rrkjus.UPF', 'O': 'O.pbe-rrkjus.UPF'}
pseudo_dir = '.'

# In this example we use a UNIX socket.  See other examples for INET socket.
# UNIX sockets are faster then INET sockets, but cannot run over a network.
# UNIX sockets are files.  The actual path will become /tmp/ipi_ase_espresso.
unixsocket = 'ase_espresso'

# Configure pw.x command for UNIX or INET.
#
# UNIX: --ipi {unixsocket}:UNIX
# INET: --ipi {host}:{port}
#
# See also QE documentation, e.g.:
#
#    https://www.quantum-espresso.org/Doc/pw_user_guide/node13.html
#
command = 'pw.x < PREFIX.pwi --ipi {unixsocket}:UNIX > PREFIX.pwo'.format(
    unixsocket=unixsocket
)

espresso = Espresso(
    command=command,
    ecutwfc=30.0,
    pseudopotentials=pseudopotentials,
    pseudo_dir=pseudo_dir,
)

opt = BFGS(atoms, trajectory='opt.traj', logfile='opt.log')

with SocketIOCalculator(
    espresso, log=sys.stdout, unixsocket=unixsocket
) as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)

# Note: QE does not generally quit cleanly - expect nonzero exit codes.
