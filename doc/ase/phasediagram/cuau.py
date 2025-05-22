# creates: cuau.png
import matplotlib.pyplot as plt

from ase.phasediagram import PhaseDiagram

refs = [
    ('Cu', 0.0),
    ('Au', 0.0),
    ('CuAu2', -0.2),
    ('CuAu', -0.5),
    ('Cu2Au', -0.7),
]
pd = PhaseDiagram(refs)
pd.plot(show=False)
plt.savefig('cuau.png')
print(pd.decompose('Cu3Au'))
