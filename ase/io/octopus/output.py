import re

import numpy as np

from ase.units import Angstrom, Bohr, Debye, Hartree, eV


class OctopusIOError(IOError):
    pass  # Cannot find output files


def read_eigenvalues_file(fd):
    unit = None

    for line in fd:
        m = re.match(r'Eigenvalues\s*\[(.+?)\]', line)
        if m is not None:
            unit = m.group(1)
            break
    line = next(fd)
    assert line.strip().startswith('#st'), line

    # fermilevel = None
    kpts = []
    eigs = []
    occs = []

    for line in fd:
        m = re.match(r'#k.*?\(\s*(.+?),\s*(.+?),\s*(.+?)\)', line)
        if m:
            k = m.group(1, 2, 3)
            kpts.append(np.array(k, float))
            eigs.append({})
            occs.append({})
        else:
            m = re.match(r'\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)', line)
            if m is None:
                m = re.match(r'Fermi energy\s*=\s*(\S+)\s*', line)
                assert m is not None
                # We can also return the fermilevel but so far we just read
                # it from the static/info instead.
                # fermilevel = float(m.group(1))
            else:
                spin, eig, occ = m.group(1, 2, 3)

                if not eigs:
                    # Only initialized if kpoint header was written
                    eigs.append({})
                    occs.append({})

                eigs[-1].setdefault(spin, []).append(float(eig))
                occs[-1].setdefault(spin, []).append(float(occ))

    nkpts = len(kpts)
    nspins = len(eigs[0])
    nbands = len(eigs[0][spin])

    kptsarr = np.array(kpts, float)
    eigsarr = np.empty((nkpts, nspins, nbands))
    occsarr = np.empty((nkpts, nspins, nbands))

    arrs = [eigsarr, occsarr]

    for arr in arrs:
        arr.fill(np.nan)

    for k in range(nkpts):
        for arr, lst in [(eigsarr, eigs), (occsarr, occs)]:
            arr[k, :, :] = [lst[k][sp] for sp
                            in (['--'] if nspins == 1 else ['up', 'dn'])]

    for arr in arrs:
        assert not np.isnan(arr).any()

    eigsarr *= {'H': Hartree, 'eV': eV}[unit]
    return kptsarr, eigsarr, occsarr


def read_static_info_stress(fd):
    stress_cv = np.empty((3, 3))

    headers = next(fd)
    assert headers.strip().startswith('T_{ij}')
    for i in range(3):
        line = next(fd)
        tokens = line.split()
        vec = np.array(tokens[1:4]).astype(float)
        stress_cv[i] = vec
    return stress_cv


def read_static_info_kpoints(fd):
    for line in fd:
        if line.startswith('List of k-points'):
            break

    tokens = next(fd).split()
    assert tokens == ['ik', 'k_x', 'k_y', 'k_z', 'Weight']
    bar = next(fd)
    assert bar.startswith('---')

    kpts = []
    weights = []

    for line in fd:
        # Format:        index   kx      ky      kz     weight
        m = re.match(r'\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)', line)
        if m is None:
            break
        kxyz = m.group(1, 2, 3)
        weight = m.group(4)
        kpts.append(kxyz)
        weights.append(weight)

    ibz_kpoints = np.array(kpts, float)
    kpoint_weights = np.array(weights, float)
    return dict(ibz_kpoints=ibz_kpoints, kpoint_weights=kpoint_weights)


def read_static_info_eigenvalues(fd, energy_unit):

    values_sknx = {}

    nbands = 0
    fermilevel = None
    for line in fd:
        line = line.strip()
        if line.startswith('#'):
            continue
        if not line[:1].isdigit():
            m = re.match(r'Fermi energy\s*=\s*(\S+)', line)
            if m is not None:
                fermilevel = float(m.group(1)) * energy_unit
            break

        tokens = line.split()
        nbands = max(nbands, int(tokens[0]))
        energy = float(tokens[2]) * energy_unit
        occupation = float(tokens[3])
        values_sknx.setdefault(tokens[1], []).append((energy, occupation))

    nspins = len(values_sknx)
    if nspins == 1:
        val = [values_sknx['--']]
    else:
        val = [values_sknx['up'], values_sknx['dn']]
    val = np.array(val, float)
    nkpts, remainder = divmod(len(val[0]), nbands)
    assert remainder == 0

    eps_skn = val[:, :, 0].reshape(nspins, nkpts, nbands)
    occ_skn = val[:, :, 1].reshape(nspins, nkpts, nbands)
    eps_skn = eps_skn.transpose(1, 0, 2).copy()
    occ_skn = occ_skn.transpose(1, 0, 2).copy()
    assert eps_skn.flags.contiguous
    d = dict(nspins=nspins,
             nkpts=nkpts,
             nbands=nbands,
             eigenvalues=eps_skn,
             occupations=occ_skn)
    if fermilevel is not None:
        d.update(fermi_level=fermilevel)
    return d


def read_static_info_energy(fd, energy_unit):
    def get(name):
        for line in fd:
            if line.strip().startswith(name):
                return float(line.split('=')[-1].strip()) * energy_unit
    return dict(energy=get('Total'), free_energy=get('Free'))


def read_static_info(fd):
    results = {}

    def get_energy_unit(line):  # Convert "title [unit]": ---> unit
        return {'[eV]': eV, '[H]': Hartree}[line.split()[1].rstrip(':')]

    for line in fd:
        if line.strip('*').strip().startswith('Brillouin zone'):
            results.update(read_static_info_kpoints(fd))
        elif line.startswith('Eigenvalues ['):
            unit = get_energy_unit(line)
            results.update(read_static_info_eigenvalues(fd, unit))
        elif line.startswith('Energy ['):
            unit = get_energy_unit(line)
            results.update(read_static_info_energy(fd, unit))
        elif line.startswith('Total stress tensor ['):
            if '[H/b^3]' in line:
                stress = read_static_info_stress(fd)
                stress *= Hartree / Bohr**3
                results.update(stress=stress)
        elif line.startswith('Total Magnetic Moment'):
            line = next(fd)
            values = line.split()
            results['magmom'] = float(values[-1])
            line = next(fd)
            assert line.startswith('Local Magnetic Moments')
            line = next(fd)
            assert line.split() == ['Ion', 'mz']
            # Reading  Local Magnetic Moments
            magmoms = []
            for line in fd:
                if line == '\n':
                    break  # there is no more thing to search for
                line = line.replace('\n', ' ')
                values = line.split()
                magmoms.append(float(values[-1]))
            results['magmoms'] = np.array(magmoms)
        elif line.startswith('Dipole'):
            assert line.split()[-1] == '[Debye]'
            dipole = [float(next(fd).split()[-1]) for i in range(3)]
            results['dipole'] = np.array(dipole) * Debye
        elif line.startswith('Forces'):
            forceunitspec = line.split()[-1]
            forceunit = {'[eV/A]': eV / Angstrom,
                         '[H/b]': Hartree / Bohr}[forceunitspec]
            forces = []
            line = next(fd)
            assert line.strip().startswith('Ion')
            for line in fd:
                if line.strip().startswith('---'):
                    break
                tokens = line.split()[-3:]
                forces.append([float(f) for f in tokens])
            results['forces'] = np.array(forces) * forceunit

    return results
