from pathlib import Path

from gpaw import GPAW, PW

from ase.build import bulk
from ase.constraints import ExpCellFilter
from ase.db import connect
from ase.dft.bandgap import bandgap
from ase.optimize import BFGS

if Path('database.db').is_file():
    Path('database.db').unlink()

structures = ['Si', 'Ge', 'C']
db = connect('database.db')

for f in structures:
    db.write(bulk(f))

for row in db.select():
    atoms = row.toatoms()
    calc = GPAW(
        mode=PW(400), kpts=(4, 4, 4), txt=f'{row.formula}-gpaw.txt', xc='LDA'
    )
    atoms.calc = calc
    atoms.get_stress()
    filter = ExpCellFilter(atoms)
    opt = BFGS(filter)
    opt.run(fmax=0.05)
    db.write(atoms=atoms, relaxed=True)


for row in db.select(relaxed=True):
    atoms = row.toatoms()
    calc = GPAW(
        mode=PW(400), kpts=(4, 4, 4), txt=f'{row.formula}-gpaw.txt', xc='LDA'
    )
    atoms.calc = calc
    atoms.get_potential_energy()
    bg, _, _ = bandgap(calc=atoms.calc)
    db.update(row.id, bandgap=bg)
