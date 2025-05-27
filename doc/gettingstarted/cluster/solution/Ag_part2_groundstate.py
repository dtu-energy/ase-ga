from gpaw import GPAW, FermiDirac

from ase.io import read

atoms = read('opt.traj')

calc = GPAW(
    mode='lcao',
    basis='sz(dzp)',
    txt='gpaw.txt',
    occupations=FermiDirac(0.1),
    setups={'Ag': '11'},
)
atoms.calc = calc
atoms.center(vacuum=4.0)
atoms.get_potential_energy()
atoms.calc.write('groundstate.gpw')
