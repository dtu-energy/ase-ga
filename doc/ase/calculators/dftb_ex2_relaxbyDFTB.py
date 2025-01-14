from ase.build import molecule
from ase.calculators.dftb import Dftb
from ase.io import read, write

atoms = molecule('H2O')
calc = Dftb(label='h2o',
            Driver_='ConjugateGradient',
            Driver_MaxForceComponent=1e-4,
            Driver_MaxSteps=1000,
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_O='p',
            Hamiltonian_MaxAngularMomentum_H='s')

atoms.calc = calc
calc.calculate(atoms)

# The 'geo_end.gen' file written by the ASE calculator
# (containing the initial geometry), has been overwritten
# by DFTB+ and now contains the final, optimized geometry.
final = read('geo_end.gen')
write('final.xyz', final)
