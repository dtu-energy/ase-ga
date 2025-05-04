from gpaw import restart

from ase.dft.wannier import Wannier

atoms, calc = restart('benzene.gpw', txt=None)

# Make wannier functions using (three) extra degrees of freedom.
wan = Wannier(nwannier=18, calc=calc, fixedstates=15)
wan.localize()
wan.save('wan18.json')
for i in range(wan.nwannier):
    wan.write_cube(i, 'benzene18_%i.cube' % i)
