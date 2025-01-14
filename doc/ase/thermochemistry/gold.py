from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.phonons import Phonons
from ase.spacegroup import crystal
from ase.thermochemistry import CrystalThermo

# Set up gold bulk and attach EMT calculator
a = 4.078
atoms = crystal('Au', (0., 0., 0.),
                spacegroup=225,
                cellpar=[a, a, a, 90, 90, 90],
                pbc=(1, 1, 1))
calc = EMT()
atoms.calc = calc
qn = QuasiNewton(atoms)
qn.run(fmax=0.05)
potentialenergy = atoms.get_potential_energy()

# Phonon analysis
N = 5
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
try:
    ph.run()
    ph.read(acoustic=True)
finally:
    ph.clean()
dos = ph.get_dos(kpts=(40, 40, 40)).sample_grid(npts=3000, width=5e-4, xmin=0.0)
phonon_energies = dos.get_energies()
phonon_DOS = dos.get_weights()


# Calculate the Helmholtz free energy
thermo = CrystalThermo(phonon_energies=phonon_energies,
                       phonon_DOS=phonon_DOS,
                       potentialenergy=potentialenergy,
                       formula_units=4)
F = thermo.get_helmholtz_energy(temperature=298.15)
