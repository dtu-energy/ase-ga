from ase.build import molecule
from ase.calculators.dftb import Dftb
from ase.io import write
from ase.optimize import QuasiNewton

atoms = molecule('H2O')
calc = Dftb(label='h2o',
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_O='p',
            Hamiltonian_MaxAngularMomentum_H='s',
            )
atoms.calc = calc

dyn = QuasiNewton(atoms, trajectory='test.traj')
dyn.run(fmax=0.01)
write('final.xyz', atoms)
