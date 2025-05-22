# creates: test.txt
import runpy
from pathlib import Path

import numpy as np

# Monkey-patch view() to avoid ASE-GUI windows popping up:
import ase.visualize

ase.visualize.view = lambda *args, **kwargs: None


def run(script):
    return runpy.run_path(Path('selfdiffusion') / script)


def run_and_get_ptp_energy(script):
    dct = run(script)
    energy = np.ptp([i.get_potential_energy() for i in dct['images']])
    return energy


e1 = run_and_get_ptp_energy('neb1.py')
assert abs(e1 - 0.111) < 0.002
e2 = run_and_get_ptp_energy('neb2.py')
assert abs(e2 - 0.564) < 0.002
e3 = run_and_get_ptp_energy('neb3.py')
assert abs(e3 - 0.239) < 0.002

with open('test.txt', 'w') as fd:
    print(e1, e2, e3, file=fd)
