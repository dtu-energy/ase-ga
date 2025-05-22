# creates: Al_phonon.png
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.phonons import Phonons

# Setup crystal and EMT calculator
atoms = bulk('Al', 'fcc', a=4.05)

# Phonon calculator
N = 7
ph = Phonons(atoms, EMT(), supercell=(N, N, N), delta=0.05)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)
ph.clean()

path = atoms.cell.bandpath('GXULGK', npoints=100)
bs = ph.get_band_structure(path)

dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)

# Plot the band structure and DOS:
import matplotlib.pyplot as plt  # noqa

fig = plt.figure(figsize=(7, 4))
ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])

emax = 0.035
bs.plot(ax=ax, emin=0.0, emax=emax)

dosax = fig.add_axes([0.8, 0.07, 0.17, 0.85])
dosax.fill_between(
    dos.get_weights(),
    dos.get_energies(),
    y2=0,
    color='grey',
    edgecolor='k',
    lw=1,
)

dosax.set_ylim(0, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel('DOS', fontsize=18)

fig.savefig('Al_phonon.png')

# --- literalinclude division line ---
# from ase.io.trajectory import Trajectory
# from ase.io import write

# Write modes for specific q-vector to trajectory files:
L = path.special_points['L']
ph.write_modes(
    [l / 2 for l in L], branches=[2], repeat=(8, 8, 8), kT=3e-4, center=True
)


# Generate gif animation:
# XXX Temporarily disabled due to matplotlib writer compatibility issue.
# with Trajectory('phonon.mode.2.traj', 'r') as traj:
#     write('Al_mode.gif', traj, interval=50,
#           rotation='-36x,26.5y,-25z')
