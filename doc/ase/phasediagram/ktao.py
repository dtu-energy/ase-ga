# creates: ktao-2d.png, ktao-3d.png
import matplotlib.pyplot as plt

from ase.phasediagram import PhaseDiagram

references = [
    ('K', 0),
    ('Ta', 0),
    ('O2', 0),
    ('K3TaO8', -16.167),
    ('KO2', -2.288),
    ('KO3', -2.239),
    ('Ta2O5', -19.801),
    ('TaO3', -8.556),
    ('TaO', -1.967),
    ('K2O', -3.076),
    ('K2O2', -4.257),
    ('KTaO3', -13.439),
]
pd = PhaseDiagram(references)
for d in [2, 3]:
    pd.plot(dims=d, show=False)
    plt.savefig(f'ktao-{d}d.png')
