from gpaw import GPAW, PW

from ase.build import bulk

atoms = bulk('Ag')
calc = GPAW(
    mode=PW(350), kpts=[8, 8, 8], txt='gpaw.bulk.Ag.txt', setups={'Ag': '11'}
)
atoms.calc = calc
atoms.get_potential_energy()
calc.write('bulk.Ag.gpw')
