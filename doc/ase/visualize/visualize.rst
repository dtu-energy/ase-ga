.. module:: ase.visualize

Visualization
=============

.. function:: view(atoms, data=None, viewer=None, repeat=None)

This provides an interface to various visualization tools, such as
:mod:`ase.gui`, RasMol_, VMD_, gOpenMol_, Avogadro_, ParaView_ or NGLView_. The default viewer is
the ase.gui, described in the :mod:`ase.gui` module. The simplest invocation
is:

.. testsetup::

    from ase import Atoms
    atoms = Atoms('Cu')

>>> from ase.visualize import view
>>> view(atoms)

where ``atoms`` is any :class:`~ase.Atoms` object.  Alternative viewers
can be used by specifying the optional keyword ``viewer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', 'rasmol', 'paraview', 'ngl'.  The VMD and Avogadro viewers can
take an optional ``data`` argument to show 3D data, such as charge density:

>>> view(atoms, viewer='VMD', data=...)

The nglview viewer additionally supports any indexible sequence of :class:`~ase.Atoms`
objects, e.g. lists of structures and :class:`~ase.io.Trajectory` objects.

If you do not wish to open an interactive gui, but rather visualize
your structure by dumping directly to a graphics file; you can use the
``write`` command of the :mod:`ase.io` module, which can write 'eps',
'png', and 'pov' files directly, like this:

>>> from ase.io import write
>>> write('image.png', atoms)

It is also possible to plot directly to a Matplotlib subplot object, which allows for
a large degree of customisability. :ref:`More information <matplotlib_plotting>`.

.. _RasMol: http://openrasmol.org/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _gOpenMol: http://www.csc.fi/gopenmol/
.. _Avogadro: http://avogadro.cc/
.. _ParaView: https://www.paraview.org/
.. _NGLView: https://github.com/arose/nglview


.. module:: ase.visualize.nglview

Viewer for Jupyter notebooks
----------------------------

A simple viewer based on X3D is built into ASE, which should work on modern
browsers without additional packages, by typing the following command into
a Jupyter_ notebook:

>>> view(atoms, viewer='x3d')

For a more feature-rich viewer, nglview is a dedicated viewer for the
Jupyter_ notebook interface.
It uses an embeddable NGL_ WebGL molecular viewer. The viewer works only in
the web browser environment and embeds the live javascript output into the
notebook. To utilize this functionality you need to have NGLView_ and
ipywidgets_ packages installed in addition to the Jupyter_ notebook.

The basic usage provided by the :func:`ase.visualize.view` function exposes
only small fraction of the NGL_ widget capabilities. The simplest form:

>>> view(atoms, viewer='ngl')

creates interactive ngl viewer widget with the few additional control widgets
added on the side. The object returned by the above call is a reference to
the `.gui` member of the :class:`ase.visualize.nglview.NGLDisplay` containing
actual viewer (`.view` member), a reference to control widgets box
(`.control_box` member) and
:func:`ase.visualize.view.nglview.NGLDisplay.custom_colors` method. The
notebook interface is not blocked by the above call and the returned object
may be further manipulated by the following code in the separate cell (the
`\_` variable contains output from the previous cell):

>>> v=_
>>> v.custom_colors({'Mn':'green','As':'blue'})
>>> v.view._remote_call("setSize", target="Widget", args=["400px", "400px"])
>>> v.view.center_view()
>>> v.view.background='#ffc'
>>> v.view.parameters=dict(clipDist=-200)

The `.view` member exposes full API of the NGLView_ widget. The
`.control_box` member is a :class:`ipywidgets.HBox` containing
:class:`nglview.widget.NGLWidget` and :class:`ipywidgets.VBox` with control
widgets. For the full documentation of these objects consult the NGLView_,
NGL_ and ipywidgets_ websites.

.. _Jupyter: https://jupyter.org/
.. _NGL: https://github.com/arose/ngl
.. _ipywidgets: https://github.com/jupyter-widgets/ipywidgets

.. autoclass:: ase.visualize.ngl.NGLDisplay
   :inherited-members:

.. autofunction:: ase.visualize.ngl.view_ngl


.. module:: ase.visualize.mlab
.. _iso surface:

Plotting iso-surfaces with Mayavi
---------------------------------

The :func:`ase.visualize.mlab.plot` function can be used from the
command-line::

    $ python -m ase.visualize.mlab abc.cube

to plot data from a cube-file or alternatively a wave function or an electron
density from a calculator restart file::

    $ python -m ase.visualize.mlab -C gpaw abc.gpw

Options:

.. include:: mlab_options.txt
    :start-after: Options:

.. autofunction:: ase.visualize.mlab.plot


.. _matplotlib_plotting:

Matplotlib
----------

>>> import matplotlib.pyplot as plt
>>> from ase.visualize.plot import plot_atoms
>>> from ase.lattice.cubic import FaceCenteredCubic
>>> slab = FaceCenteredCubic('Au', size=(2, 2, 2))
>>> fig, ax = plt.subplots()
>>> plot_atoms(slab, ax, radii=0.3, rotation=('90x,45y,0z'))
>>> fig.savefig("ase_slab.png")

.. image:: matplotlib_plot_atoms1.png

The data is plotted directly to a Matplotlib subplot object, giving a large
degree of customizability.


>>> import matplotlib.pyplot as plt
>>> from ase.visualize.plot import plot_atoms
>>> from ase.lattice.cubic import FaceCenteredCubic
>>> slab = FaceCenteredCubic('Au', size=(2, 2, 2))
>>> fig, axarr = plt.subplots(1, 4, figsize=(15, 5))
>>> plot_atoms(slab, axarr[0], radii=0.3, rotation=('0x,0y,0z'))
>>> plot_atoms(slab, axarr[1], scale=0.7, offset=(3, 4), radii=0.3, rotation=('0x,0y,0z'))
>>> plot_atoms(slab, axarr[2], radii=0.3, rotation=('45x,45y,0z'))
>>> plot_atoms(slab, axarr[3], radii=0.3, rotation=('0x,0y,0z'))
>>> axarr[0].set_title("No rotation")
>>> axarr[1].set_xlabel("X-axis, [$\mathrm{\AA}$]")
>>> axarr[1].set_ylabel("Y-axis, [$\mathrm{\AA}$]")
>>> axarr[2].set_axis_off()
>>> axarr[3].set_xlim(2, 6)
>>> axarr[3].set_ylim(2, 6)
>>> fig.savefig("ase_slab_multiple.png")

.. image:: matplotlib_plot_atoms2.png

>>> stem_image = mpimg.imread("stem_image.jpg")
>>> atom_pos = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5), (0.5, 0.5, 0.0)]
>>> srtio3 = crystal(['Sr','Ti','O'], atom_pos, spacegroup=221, cellpar=3.905, size=(3, 3, 3))
>>> fig, ax = plt.subplots()
>>> ax.imshow(stem_image, cmap='gray')
>>> plot_atoms(srtio3, ax, radii=0.3, scale=6.3, offset=(47, 54), rotation=('90x,45y,56z'))
>>> ax.set_xlim(0, stem_image.shape[0])
>>> ax.set_ylim(0, stem_image.shape[1])
>>> ax.set_axis_off()
>>> fig.tight_layout()
>>> fig.savefig("iomatplotlib3.png")

.. image:: matplotlib_plot_atoms3.png
