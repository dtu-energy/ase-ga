from ase.db import connect

refs = connect('refs.db')
db = connect('ads.db')

for row in db.select():
    ea = (
        row.energy
        - refs.get(formula=row.ads).energy
        - refs.get(layers=row.layers, surf=row.surf).energy
    )
    h = row.positions[-1, 2] - row.positions[-2, 2]
    db.update(row.id, height=h, ea=ea)
