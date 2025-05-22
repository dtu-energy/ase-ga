from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations

atoms = molecule('N2')
atoms.calc = EMT()
dyn = QuasiNewton(atoms)
dyn.run(fmax=0.01)
potentialenergy = atoms.get_potential_energy()

vib = Vibrations(atoms)
vib.run()
vib_energies = vib.get_energies()

thermo = IdealGasThermo(
    vib_energies=vib_energies,
    potentialenergy=potentialenergy,
    atoms=atoms,
    geometry='linear',
    symmetrynumber=2,
    spin=0,
)
G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.0)
