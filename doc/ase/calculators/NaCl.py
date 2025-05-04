from ase import Atom, Atoms
from ase.calculators.vasp import Vasp

a = [6.5, 6.5, 7.7]
d = 2.3608
NaCl = Atoms(
    [Atom('Na', [0, 0, 0], magmom=1.928), Atom('Cl', [0, 0, d], magmom=0.75)],
    cell=a,
)

calc = Vasp(prec='Accurate', xc='PBE', lreal=False)
NaCl.calc = calc

print(NaCl.get_magnetic_moment())
