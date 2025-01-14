# creates: score-size-fcc2sc.svg

import glob
import json
import os

import matplotlib.pyplot as plt

for fname in glob.glob('Popt*2*.json'):
    tag = os.path.basename(fname).replace('Popt-', '').replace('.json', '')

    with open(fname) as data_file:
        data = json.load(data_file)
    x = []
    y = []
    for nuc, rec in sorted(data.items()):
        x.append(int(nuc))
        y.append(rec['dev'])

    plt.figure(figsize=(4, 3))
    plt.text(1950, 0.35,
             tag.replace('2', r' $\rightarrow$ '), horizontalalignment='right')
    plt.xlabel(r'Number of primitive unit cells $N_{uc}$')
    plt.ylabel(r'Optimality measure $\bar \Delta$')
    plt.axis([0, 2000, -0.025, 0.4])
    plt.plot(x, y, 'bo')
    plt.savefig('score-size-%s.svg' % tag, bbox_inches='tight')
