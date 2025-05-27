from ase.build import bulk
from ase.calculators.emt import EMT
from ase.db import connect
from ase.eos import calculate_eos

db = connect('bulk.db')
for symb in ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']:
    atoms = bulk(symb, 'fcc')
    atoms.calc = EMT()
    eos = calculate_eos(atoms)
    v, e, B = eos.fit()  # find minimum
    # Do one more calculation at the minimu and write to database:
    atoms.cell *= (v / atoms.get_volume()) ** (1 / 3)
    atoms.get_potential_energy()
    db.write(atoms, bm=B)
