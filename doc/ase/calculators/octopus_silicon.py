from ase.build import bulk
from ase.calculators.octopus import Octopus

system = bulk('Si')

calc = Octopus(
    directory='oct-si',
    Spacing=0.25,
    KPointsGrid=[[4, 4, 4]],
    KPointsUseSymmetries=True,
    Output=[['dos'], ['density'], ['potential']],
    OutputFormat='xcrysden',
    DosGamma=0.1,
)

system.calc = calc
system.get_potential_energy()
