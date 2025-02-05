import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(dpi=80, figsize=(4.2, 6))
fig.subplots_adjust(left=.16, right=.97, top=.97, bottom=.05)

# Plot KS bands
k, eps = np.loadtxt('KSbands.txt', unpack=True)
plt.plot(k, eps, 'ro', label='DFT', ms=9)

# Plot Wannier bands
k, eps = np.loadtxt('WANbands.txt', unpack=True)
plt.plot(k, eps, 'k.', label='Wannier')

plt.plot([-.5, .5], [1, 1], 'k:', label='_nolegend_')
plt.text(-.5, 1, 'fixedenergy', ha='left', va='bottom')
plt.axis('tight')
plt.xticks([-.5, -.25, 0, .25, .5],
           [r'$X$', r'$\Delta$', r'$\Gamma$', r'$\Delta$', r'$X$'], size=16)
plt.ylabel(r'$E - E_F\  \rm{(eV)}$', size=16)
plt.legend()
plt.savefig('bands.png', dpi=80)
plt.show()
