from gpaw import restart

from ase.dft.wannier import Wannier

atoms, calc = restart('benzene.gpw', txt=None)

# Make wannier functions of occupied space only
wan = Wannier(nwannier=15, calc=calc)
wan.localize()
for i in range(wan.nwannier):
    wan.write_cube(i, 'benzene15_%i.cube' % i)
