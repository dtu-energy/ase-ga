# fmt: off
import numpy as np
import pytest

from ase.build import bulk
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.io import read, write

magres_with_complex_labels = """#$magres-abinitio-v1.0
# Fictional test data
[atoms]
units lattice Angstrom
lattice       8.0 0.0 0.0 0.0 10.2 0.0 0.0 0.0 15.9
units atom Angstrom
atom H H1   1    3.70 4.9 1.07
atom H H2   100  4.30 2.1 5.15
atom H H2a  101  4.30 2.1 5.15
atom H H2b  102  4.30 2.1 5.15
[/atoms]
[magres]
units ms ppm
ms H1  1    3.0 -5.08 -3.19 -4.34 3.42 -3.78  1.05 -1.89  2.85
ms H2100    3.0  5.08  3.19  4.34 3.42 -3.78 -1.05 -1.89  2.85
ms H2a101   3.0  5.08  3.19  4.34 3.42 -3.78 -1.05 -1.89  2.85
ms H2b 102  3.0  5.08  3.19  4.34 3.42 -3.78 -1.05 -1.89  2.85
units efg au
efg H1  1   9.6 -3.43  1.45 -3.43 8.52 -1.43  1.45 -1.43 -1.81
efg H2100   9.6  3.43 -1.45  3.43 8.52 -1.43 -1.45 -1.43 -1.81
efg H2a101  9.6  3.43 -1.45  3.43 8.52 -1.43 -1.45 -1.43 -1.81
efg H2b 102 9.6  3.43 -1.45  3.43 8.52 -1.43 -1.45 -1.43 -1.81
[/magres]
"""

magres_with_too_large_index = """#$magres-abinitio-v1.0
# Test data with index >999
[atoms]
units lattice Angstrom
lattice       10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0
units atom Angstrom
atom H H1 1000    0.0 0.0 0.0
[/atoms]
[magres]
units ms ppm
ms H11000    1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
[/magres]
"""


def test_magres():

    # Test with fictional data
    si2 = bulk('Si')

    ms = np.ones((2, 3, 3))
    si2.set_array('ms', ms)
    efg = np.repeat([[[1, 0, 0], [0, 1, 0], [0, 0, -2]]], 2, axis=0)
    si2.set_array('efg', efg)

    calc = SinglePointDFTCalculator(si2)
    calc.results['sus'] = np.eye(3) * 2
    si2.calc = calc

    si2.info['magres_units'] = {'ms': 'ppm',
                                'efg': 'au',
                                'sus': '10^-6.cm^3.mol^-1'}

    write('si2_test.magres', si2)
    si2 = read('si2_test.magres')

    assert (np.trace(si2.get_array('ms')[0]) == 3)
    assert (np.all(np.isclose(si2.get_array('efg')[:, 2, 2], -2)))
    assert (np.all(np.isclose(si2.calc.results['sus'], np.eye(3) * 2)))


def test_magres_large(datadir):

    # Test with big structure
    assert len(read(datadir / "large_atoms.magres")) == 240


def test_magres_sitelabels(datadir):
    """Test reading magres files with munged site labels and indices
    for cases where the site label has a number in it."""

    # Write temporary file
    with open('magres_with_complex_labels.magres', 'w') as f:
        f.write(magres_with_complex_labels)

    # Read it back
    atoms = read('magres_with_complex_labels.magres')

    labels_ref = ['H1', 'H2', 'H2a', 'H2b']
    labels = atoms.get_array('labels')
    np.testing.assert_array_equal(labels, labels_ref)

    indices_ref = [1, 100, 101, 102]
    indices = atoms.get_array('indices')
    np.testing.assert_array_equal(indices, indices_ref)


def test_magres_with_large_indices():
    """Test handling of magres files with indices >999"""
    # Write temporary file
    with open('magres_large_index.magres', 'w') as f:
        f.write(magres_with_too_large_index)

    # Check that reading raises the correct error
    with pytest.raises(RuntimeError, match="Index greater than 999 detected"):
        read('magres_large_index.magres')
