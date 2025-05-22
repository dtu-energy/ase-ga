from ase.io.opls import OPLSStructure
from ase.visualize import view

# 172_mod.xyz if the file name for the structure above
s = OPLSStructure('172_mod.xyz')
view(s)  # view with real elements
elements = {'CT': 'Si', 'HC': 'H', 'H1': 'He'}
view(s.colored(elements))  # view with fake elements
