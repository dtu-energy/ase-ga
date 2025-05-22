from gpaw import GPAW

calc = GPAW('bulk.Ag.gpw')
atoms = calc.get_atoms()
path = atoms.cell.bandpath('WLGXWK', density=10)
path.write('path.json')

calc = calc.fixed_density(kpts=path, symmetry='off')

bs = calc.band_structure()
bs.write('bs.json')
