# creates: ptable.png
from ase.utils.ptable import ptable

atoms = ptable()
atoms.write('ptable.png')

#   from ase.visualize import view
#   view(atoms)
