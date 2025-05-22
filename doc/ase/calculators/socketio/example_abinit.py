from ase.build import bulk
from ase.calculators.abinit import Abinit
from ase.calculators.socketio import SocketIOCalculator
from ase.constraints import ExpCellFilter
from ase.optimize import BFGS

atoms = bulk('Si')
atoms.rattle(stdev=0.1, seed=42)

# Configuration parameters; please edit as appropriate
pps = '/path/to/pseudopotentials'
pseudopotentials = {'Si': '14-Si.LDA.fhi'}
exe = 'abinit'

unixsocket = 'ase_abinit'
command = f'{exe} PREFIX.in --ipi {unixsocket}:UNIX > PREFIX.log'
# (In the command, note that PREFIX.in must precede --ipi.)


configuration_kwargs = dict(
    command=command,
    pp_paths=[pps],
)


# Implementation note: Socket-driven calculations in Abinit inherit several
# controls for from the ordinary cell optimization code.  We have to hack those
# variables in order for Abinit not to decide that the calculation converged:
boilerplate_kwargs = dict(
    ionmov=28,  # activate i-pi/socket mode
    expert_user=1,  # Ignore warnings (chksymbreak, chksymtnons, chkdilatmx)
    optcell=2,  # allow the cell to relax
    tolmxf=1e-300,  # Prevent Abinit from thinking we "converged"
    ntime=100_000,  # Allow "infinitely" many iterations in Abinit
    ecutsm=0.5,  # Smoothing PW cutoff energy (mandatory for cell optimization)
)


kwargs = dict(
    ecut=5 * 27.3,
    tolvrs=1e-8,
    kpts=[2, 2, 2],
    **boilerplate_kwargs,
    **configuration_kwargs,
)

abinit = Abinit(**kwargs)

opt = BFGS(ExpCellFilter(atoms), trajectory='opt.traj')

with SocketIOCalculator(abinit, unixsocket=unixsocket) as atoms.calc:
    opt.run(fmax=0.01)
