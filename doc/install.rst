.. _download_and_install:

============
Installation
============

Requirements
============

* Python_ 3.9 or newer
* NumPy_ (base N-dimensional array package)
* SciPy_ (library for scientific computing)
* Matplotlib_ (plotting)

Optional:

* Flask_ for :mod:`ase.db` web-interface
* pytest_ for running tests
* pytest-mock_ for running some more tests
* pytest-xdist_ for running tests in parallel
* spglib_ for certain symmetry-related features

.. _Python: https://www.python.org/
.. _NumPy: https://docs.scipy.org/doc/numpy/reference/
.. _SciPy: https://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: https://matplotlib.org/
.. _Flask: https://palletsprojects.com/p/flask/
.. _PyPI: https://pypi.org/project/ase
.. _PIP: https://pip.pypa.io/en/stable/
.. _pytest: https://pypi.org/project/pytest/
.. _pytest-xdist: https://pypi.org/project/pytest-xdist/
.. _pytest-mock: https://pypi.org/project/pytest-mock/
.. _spglib: https://pypi.org/project/spglib/

Installation using system package managers
==========================================

Linux
-----

Major GNU/Linux distributions (including Debian and Ubuntu derivatives,
Arch, Fedora, Red Hat and CentOS) have a ``python-ase`` package
available that you can install on your system. This will manage
dependencies and make ASE available for all users.

.. note::
   Depending on the distribution, this may not be the latest
   release of ASE.

Max OSX (Homebrew)
------------------

It is generally not recommended to rely on the Mac system Python.
Mac OSX versions before 12.3 included an old Python, incompatible with
ASE and missing the pip_ package manager. Since 12.3 MacOS has moved
to newer Python versions but it may not be installed by default.
Over time, different approaches to Python on Mac have been popular and things can get a bit messy_.

Before installing ASE with ``pip`` as described in the next section, Mac
users need to install an appropriate Python version.
One approach uses the Homebrew_ package manager, which provides an up-to-date version
of Python 3 and the tkinter library needed for ``ase gui``::

  $ brew install python-tk

Note that Homebrew only allows ``pip`` to install to virtual environments.
For more information about the quirks of brewed Python see this guide_.

.. _messy: https://xkcd.com/1987/

.. _Homebrew: http://brew.sh

.. _guide: https://docs.brew.sh/Homebrew-and-Python


.. index:: pip
.. _pip installation:


Installation from PyPI using pip
================================

.. highlight:: bash

The simplest way to install ASE is to use pip_ which will automatically get
the source code from PyPI_::

    $ pip install --upgrade ase

If you intend to `run the tests`_, use::

    $ pip install --upgrade ase[test]


Local user installation
-----------------------

The above commands should work if you are using a virtualenv, Conda or
some other user-level Python installation.

If your Python/pip was provided by a system administrator or package
manager, you might not have appropriate permissions to install ASE and
its dependencies to the default locations. In that case install with e.g.::

    $ pip install --upgrade --user ase

This will install ASE in a local folder where Python can
automatically find it (``~/.local`` on Unix, see here_ for details).  The
:ref:`cli` will be installed in the following location:

=================  ============================
Unix and Mac OS X  ``~/.local/bin``
Windows            ``%APPDATA%/Python/Scripts``
=================  ============================

Make sure you have that path in your :envvar:`PATH` environment variable.

.. _here: https://docs.python.org/3/library/site.html#site.USER_BASE

.. _download:

Installation from source
========================

As an alternative to PyPI, you can also get the source from a tar-file or
from Git.

:Tar-file:

    You can get the source as a `tar-file <http://xkcd.com/1168/>`__ for the
    latest stable release (ase-3.23.0.tar.gz_) or the latest
    development snapshot (`<snapshot.tar.gz>`_).

    Unpack and make a soft link::

        $ tar -xf ase-3.23.0.tar.gz
        $ ln -s ase-3.23.0 ase

    Here is a `list of tarballs <https://pypi.org/simple/ase/>`__.

:Git clone:

    Alternatively, you can get the source for the latest stable release from
    https://gitlab.com/ase/ase like this::

        $ git clone -b 3.23.0 https://gitlab.com/ase/ase.git

    or if you want the development version::

        $ git clone https://gitlab.com/ase/ase.git


With the source from a Git clone or tar file, you can install the code with ``pip install /path/to/source``,
which will manage dependencies as though installing from PyPI.
(See `Local user installation`_ above if there are permissions problems.)
Alternatively, you can add ``~/ase`` to your :envvar:`PYTHONPATH` environment variable
and add ``~/ase/bin`` to :envvar:`PATH` (assuming ``~/ase`` is where your ASE folder is).
In this case you are responsible for also installing the dependencies listed in `pyproject.toml`_.

Finally, please `run the tests`_.

.. _pyproject.toml : https://gitlab.com/ase/ase/-/blob/master/pyproject.toml


Pip install directly from git source
------------------------------------

This is a convenient way to install the "bleeding-edge" master
branch directly with pip, if you don't intend to do further development::

    $ pip install --upgrade git+https://gitlab.com/ase/ase.git@master

The ``--upgrade`` ensures that you always reinstall even if the version
number hasn't changed.


.. note::

    We also have Git-tags for older stable versions of ASE.
    See the :ref:`releasenotes` for which tags are available.  Also the
    dates of older releases can be found there.


.. _ase-3.23.0.tar.gz: https://pypi.org/packages/source/a/ase/ase-3.23.0.tar.gz

.. index:: test
.. _running tests:
.. _run the tests:

Test your installation
======================

Before running the tests, make sure you have set your :envvar:`PATH`
environment variable correctly as described in the relevant section above.
Run the tests like this::

    $ ase test  # takes 1 min.

and send us the output if there are failing tests.
