from ase.spacegroup import crystal

a = 3.21
c = 5.21
mg = crystal(
    'Mg',
    [(1.0 / 3.0, 2.0 / 3.0, 3.0 / 4.0)],
    spacegroup=194,
    cellpar=[a, a, c, 90, 90, 120],
)
