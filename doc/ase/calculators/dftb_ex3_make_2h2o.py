# fun collision of:  2 H2 + O2 -> 2 H2O
import os

from ase.build import molecule
from ase.calculators.dftb import Dftb
from ase.io import read, write
from ase.io.dftb import read_dftb_velocities, write_dftb_velocities

o2 = molecule('O2')
h2_1 = molecule('H2')
h2_2 = molecule('H2')
o2.translate([0, 0.01, 0])
h2_1.translate([0, 0, 3])
h2_1.euler_rotate(center='COP', theta=90)
h2_2.translate([0, 0, -3])
h2_2.euler_rotate(center='COP', theta=90)
o2.set_velocities(([0, 0, 0], [0, 0, 0]))
h2_1.set_velocities(([0, 0, -3.00], [0, 0, -3.000]))
h2_2.set_velocities(([0, 0, 3.000], [0, 0, 3.000]))
atoms = o2 + h2_1 + h2_2

# 1fs = 41.3 au
# 1000K = 0.0031668 au
calculator_NVE = Dftb(label='h2o',
                      Hamiltonian_MaxAngularMomentum_='',
                      Hamiltonian_MaxAngularMomentum_O='p',
                      Hamiltonian_MaxAngularMomentum_H='s',
                      Driver_='VelocityVerlet',
                      Driver_MDRestartFrequency=10,
                      Driver_Velocities_='',
                      Driver_Velocities_empty='<<+ "velocities.txt"',
                      Driver_Steps=1000,
                      Driver_KeepStationary='Yes',
                      Driver_TimeStep=4.13,
                      Driver_Thermostat_='None',
                      Driver_Thermostat_empty='')

# 1fs = 41.3 au
# 1000K = 0.0031668 au
calculator_NVT = Dftb(label='h2o',
                      Hamiltonian_MaxAngularMomentum_='',
                      Hamiltonian_MaxAngularMomentum_O='p',
                      Hamiltonian_MaxAngularMomentum_H='s',
                      Driver_='VelocityVerlet',
                      Driver_MDRestartFrequency=5,
                      Driver_Velocities_='',
                      Driver_Velocities_empty='<<+ "velocities.txt"',
                      Driver_Steps=500,
                      Driver_KeepStationary='Yes',
                      Driver_TimeStep=8.26,
                      Driver_Thermostat_='Berendsen',
                      Driver_Thermostat_Temperature=0.00339845142,  # 800 degC
                      Driver_Thermostat_CouplingStrength=0.01)

write_dftb_velocities(atoms, 'velocities.txt')

atoms.calc = calculator_NVE
atoms.get_potential_energy()  # run NVE ensemble using DFTB+'s own driver
atoms = read('geo_end.gen')
write('after_NVE.xyz', atoms)

read_dftb_velocities(atoms, filename='geo_end.xyz')
write_dftb_velocities(atoms, 'velocities.txt')
os.system('mv geo_end.xyz geo_end_NVE.xyz')

atoms.calc = calculator_NVT
atoms.get_potential_energy()  # run NVT ensemble using DFTB+'s own driver
atoms = read('geo_end.gen')
write('after_NVT.xyz', atoms)

read_dftb_velocities(atoms, filename='geo_end.xyz')
write_dftb_velocities(atoms, 'velocities.txt')
os.system('mv geo_end.xyz geo_end_NVT.xyz')

# to watch:
#  ase gui geo_end_NVE.xyz geo_end_NVT.xyz
