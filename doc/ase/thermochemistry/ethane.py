from numpy import array

from ase.thermochemistry import HinderedThermo

vibs = array(
    [
        3049.060670,
        3040.796863,
        3001.661338,
        2997.961647,
        2866.153162,
        2750.855460,
        1436.792655,
        1431.413595,
        1415.952186,
        1395.726300,
        1358.412432,
        1335.922737,
        1167.009954,
        1142.126116,
        1013.918680,
        803.400098,
        783.026031,
        310.448278,
        136.112935,
        112.939853,
        103.926392,
        77.262869,
        60.278004,
        25.825447,
    ]
)
vib_energies = vibs / 8065.54429  # convert to eV from cm^-1
trans_barrier_energy = 0.049313  # eV
rot_barrier_energy = 0.017675  # eV
sitedensity = 1.5e15  # cm^-2
rotationalminima = 6
symmetrynumber = 1
mass = 30.07  # amu
inertia = 73.149  # amu Ang^-2

thermo = HinderedThermo(
    vib_energies=vib_energies,
    trans_barrier_energy=trans_barrier_energy,
    rot_barrier_energy=rot_barrier_energy,
    sitedensity=sitedensity,
    rotationalminima=rotationalminima,
    symmetrynumber=symmetrynumber,
    mass=mass,
    inertia=inertia,
)

F = thermo.get_helmholtz_energy(temperature=298.15)
