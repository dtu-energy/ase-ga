from gpaw import GPAW

calc = GPAW('groundstate.rutile.gpw')
atoms = calc.get_atoms()
path = atoms.cell.bandpath(density=7)
path.write('path.rutile.json')

calc = calc.fixed_density(kpts=path, symmetry='off')

bs = calc.band_structure()
bs.write('bs.rutile.json')
