from ase.build import niggli_reduce
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import ExpCellFilter
from ase.ga import set_raw_score
from ase.optimize import FIRE

try:
    from asap3 import EMT
except ImportError:
    from ase.calculators.emt import EMT


def finalize(atoms, energy=None, forces=None, stress=None):
    # Finalizes the atoms by attaching a SinglePointCalculator
    # and setting the raw score as the negative of the total energy
    atoms.wrap()
    calc = SinglePointCalculator(
        atoms, energy=energy, forces=forces, stress=stress
    )
    atoms.calc = calc
    raw_score = -atoms.get_potential_energy()
    set_raw_score(atoms, raw_score)


def relax(atoms, cellbounds=None):
    # Performs a variable-cell relaxation of the structure
    calc = EMT()
    atoms.calc = calc

    converged = False
    niter = 0
    while not converged and niter < 10:
        if cellbounds is not None:
            cell = atoms.get_cell()
            if not cellbounds.is_within_bounds(cell):
                niggli_reduce(atoms)
            cell = atoms.get_cell()
            if not cellbounds.is_within_bounds(cell):
                # Niggli reduction did not bring the unit cell
                # within the specified bounds; this candidate should
                # be discarded so we set an absurdly high energy
                finalize(atoms, 1e9)
                return

        ecf = ExpCellFilter(atoms)
        dyn = FIRE(ecf, maxmove=0.2, logfile=None, trajectory=None)
        dyn.run(fmax=1e-3, steps=100)

        converged = dyn.converged()
        niter += 1

    dyn = FIRE(atoms, maxmove=0.2, logfile=None, trajectory=None)
    dyn.run(fmax=1e-2, steps=100)

    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    s = atoms.get_stress()
    finalize(atoms, energy=e, forces=f, stress=s)
