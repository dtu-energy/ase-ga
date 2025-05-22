from gpaw import GPAW

from ase import Atoms

atoms = Atoms('N')
atoms.center(vacuum=3.0)
atoms.set_initial_magnetic_moments([3])

calc = GPAW(mode='lcao', basis='dzp')
atoms.calc = calc
atoms.get_potential_energy()
