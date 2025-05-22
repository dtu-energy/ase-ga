import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(dpi=80, figsize=(4.2, 6))
fig.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.05)

# Plot KS bands
k, eps = np.loadtxt('KSbands.txt', unpack=True)
plt.plot(k, eps, 'ro', label='DFT', ms=9)

# Plot Wannier bands
k, eps = np.loadtxt('WANbands.txt', unpack=True)
plt.plot(k, eps, 'k.', label='Wannier')

plt.plot([-0.5, 0.5], [1, 1], 'k:', label='_nolegend_')
plt.text(-0.5, 1, 'fixedenergy', ha='left', va='bottom')
plt.axis('tight')
plt.xticks(
    [-0.5, -0.25, 0, 0.25, 0.5],
    [r'$X$', r'$\Delta$', r'$\Gamma$', r'$\Delta$', r'$X$'],
    size=16,
)
plt.ylabel(r'$E - E_F\  \rm{(eV)}$', size=16)
plt.legend()
plt.savefig('bands.png', dpi=80)
plt.show()
