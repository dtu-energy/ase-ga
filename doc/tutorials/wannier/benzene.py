from gpaw import GPAW

from ase.build import molecule

atoms = molecule('C6H6')
atoms.center(vacuum=3.5)

calc = GPAW(mode='fd', h=0.21, xc='PBE', txt='benzene.txt', nbands=18)
atoms.calc = calc
atoms.get_potential_energy()

calc = calc.fixed_density(
    txt='benzene-harris.txt',
    nbands=40,
    eigensolver='cg',
    convergence={'bands': 35},
)
atoms.get_potential_energy()

calc.write('benzene.gpw', mode='all')
