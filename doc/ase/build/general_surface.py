# creates: s1.png s2.png s3.png s4.png general_surface.pdf
import os
import shutil
from pathlib import Path

from ase import Atoms
from ase.build import bulk, surface
from ase.io import write

# Au example:
s1 = surface('Au', (2, 1, 1), 9)
s1.center(vacuum=10, axis=2)
# Mo example:
Mobulk = bulk('Mo', 'bcc', a=3.16, cubic=True)
s2 = surface(Mobulk, (3, 2, 1), 9)
s2.center(vacuum=10, axis=2)
# Pt3Rh example:
a = 4.0
Pt3Rh = Atoms('Pt3Rh',
              scaled_positions=[(0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5)],
              cell=[a, a, a],
              pbc=True)
s3 = surface(Pt3Rh, (2, 1, 1), 9)
s3.center(vacuum=10, axis=2)

Pt3Rh.set_chemical_symbols('PtRhPt2')
s4 = surface(Pt3Rh, (2, 1, 1), 9)
s4.center(vacuum=10, axis=2)
# Done

for atoms, name in [(s1, 's1'), (s2, 's2'), (s3, 's3'), (s4, 's4')]:
    write(name + '.pov', atoms,
          rotation='-90x',
          povray_settings=dict(
              transparent=False)).render()


dir = os.environ.get('PDF_FILE_DIR')
if dir:
    shutil.copyfile(Path(dir) / 'general_surface.pdf',
                    'general_surface.pdf')
else:
    for i in range(2):
        error = os.system(
            'pdflatex -interaction=nonstopmode general_surface > /dev/null')
        if error:
            with open('general_surface.pdf', 'w') as fd:
                fd.write('pdflatex not found\n')
            break
        os.remove('general_surface.aux')
        os.remove('general_surface.log')
