Atomic Simulation Environment
=============================

ASE is a set of tools and Python modules for setting up, manipulating,
running, visualizing and analyzing atomistic simulations.

Webpage: http://wiki.fysik.dtu.dk/ase


Requirements
------------

* Python_ 3.9 or later
* NumPy_ (base N-dimensional array package)
* SciPy_ (library for scientific computing)
* Matplotlib_ (2D Plotting)

Optional:

* Flask_ (for ase.db web-interface)
* spglib_ (for symmetry operations)

Installation
------------

To install the latest release:

::

  pip install ase

To install the development version:

::

  pip install git+https://gitlab.com/ase/ase.git

Testing
-------

Please run the tests::

    $ ase test  # takes 1 min.

and send us the output if there are failing tests.


Contact
-------

* Mailing list: ase-users_

* Chat: Join the ``#ase`` channel on Matrix_, also accessible via the Element_ webclient.

* There is an `ASE forum <https://matsci.org/c/ase/36>`_ on
  the `Materials Science Community Forum <https://matsci.org/>`_.

Feel free to create Merge Requests and Issues on our GitLab page:
https://gitlab.com/ase/ase

For regular support, please use the mailing list or chat rather than GitLab.


Example
-------

Geometry optimization of hydrogen molecule with NWChem:

>>> from ase import Atoms
>>> from ase.optimize import BFGS
>>> from ase.calculators.nwchem import NWChem
>>> from ase.io import write
>>> h2 = Atoms('H2',
               positions=[[0, 0, 0],
                          [0, 0, 0.7]])
>>> h2.calc = NWChem(xc='PBE')
>>> opt = BFGS(h2, trajectory='h2.traj')
>>> opt.run(fmax=0.02)
BFGS:   0  19:10:49    -31.435229     2.2691
BFGS:   1  19:10:50    -31.490773     0.3740
BFGS:   2  19:10:50    -31.492791     0.0630
BFGS:   3  19:10:51    -31.492848     0.0023
>>> write('H2.xyz', h2)
>>> h2.get_potential_energy()  # ASE's units are eV and Ang
-31.492847800329216

This example requires NWChem to be installed.

::

    $ ase gui h2.traj


.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _flask: https://pypi.org/project/Flask/
.. _spglib: https://github.com/spglib/spglib
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _Matrix: https://matrix.to/#/!JEiuNJLuxedbohAOuH:matrix.org
.. _Element: https://app.element.io/#/room/#ase:matrix.org
