# creates: formats.csv
from ase.formula import Formula

formulas = ['H2O', 'SiC', 'MoS2', 'AB2', 'BN', 'SiO2']
formats = [
    'hill',
    'metal',
    'abc',
    'reduce',
    'ab2',
    'a2b',
    'periodic',
    'latex',
    'html',
    'rest',
]


with open('formats.csv', 'w') as fd:
    print(', '.join(formats), file=fd)
    for f in formulas:
        formula = Formula(f)
        print(', '.join(formula.format(fmt) for fmt in formats), file=fd)
