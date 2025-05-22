from ase import Atoms
from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.db import connect

db1 = connect('bulk.db')
db2 = connect('ads.db')


def run(symb, a, n):
    atoms = fcc111(symb, (1, 1, n), a=a)
    atoms.calc = EMT()
    atoms.get_forces()
    return atoms


# Clean slabs:
for row in db1.select():
    a = row.cell[0, 1] * 2
    symb = row.symbols[0]
    for n in [1, 2, 3]:
        id = db2.reserve(layers=n, surf=symb, ads='clean')
        if id is not None:
            atoms = run(symb, a, n)
            db2.write(atoms, id=id, layers=n, surf=symb, ads='clean')

# Atoms:
for ads in 'CNO':
    a = Atoms(ads)
    a.calc = EMT()
    a.get_potential_energy()
    db2.write(a)
