import numpy as np
from gpaw import restart

from ase.dft.wannier import Wannier

atoms, calc = restart('poly.gpw', txt=None)

# Make wannier functions using (one) extra degree of freedom
wan = Wannier(
    nwannier=6,
    calc=calc,
    fixedenergy=1.5,
    initialwannier='orbitals',
    functional='var',
)
wan.localize()
wan.save('poly.json')
wan.translate_all_to_cell((2, 0, 0))
for i in range(wan.nwannier):
    wan.write_cube(i, 'polyacetylene_%i.cube' % i)

# Print Kohn-Sham bandstructure
ef = calc.get_fermi_level()
with open('KSbands.txt', 'w') as fd:
    for k, kpt_c in enumerate(calc.get_ibz_k_points()):
        for eps in calc.get_eigenvalues(kpt=k):
            print(kpt_c[0], eps - ef, file=fd)

# Print Wannier bandstructure
with open('WANbands.txt', 'w') as fd:
    for k in np.linspace(-0.5, 0.5, 100):
        ham = wan.get_hamiltonian_kpoint([k, 0, 0])
        for eps in np.linalg.eigvalsh(ham).real:
            print(k, eps - ef, file=fd)
