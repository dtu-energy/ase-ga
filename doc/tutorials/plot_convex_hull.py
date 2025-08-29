import matplotlib.pyplot as plt
import numpy as np

from ase.db import connect
from ase.io import write
from ase.phasediagram import PhaseDiagram

db = connect('hull.db')

# Select the evaluated candidates and retrieve the chemical formula and mixing
# energy for the phase diagram
refs = []
dcts = list(db.select('relaxed=1'))
for dct in dcts:
    refs.append((dct.formula, -dct.raw_score))

pd = PhaseDiagram(refs)
ax = pd.plot(
    show=not True,  # set to True to show plot
    only_label_simplices=True,
)
plt.savefig('hull.png')

# View the simplices of the convex hull
simplices = []
toview = sorted(np.array(dcts)[pd.hull], key=lambda x: x.mass)
for dct in toview:
    simplices.append(dct.toatoms())

write('hull.traj', simplices)
