.. module:: ase.calculators.kim

===
KIM
===

.. note::
   This package requires the *KIM API* package, which is `hosted on GitHub
   <https://github.com/openkim/kim-api>`__ and available through many binary package
   managers.  See `openkim.org/kim-api <https://openkim.org/kim-api>`_ for installation
   options.

.. note::
   This package requires the *kimpy* python package, which is `hosted on GitHub
   <https://github.com/openkim/kimpy>`__ and also made available through `PyPI
   <https://pypi.org/project/kimpy/>`_.

.. _Overview:

--------
Overview
--------

This package contains a calculator interface that allows one to easily use any potential
archived in `Open Knowledgebase of Interatomic Models (OpenKIM) <https://openkim.org>`_
through ASE.  OpenKIM is an NSF-funded project aimed at providing easy access to
standardized implementations of classical interatomic potentials that can be used with a
variety of molecular simulation codes.

If you haven't done so already, you'll need to install the `KIM Application Programming
Interface (API) <https://openkim.org/kim-api>`_ and the `kimpy
<https://pypi.org/project/kimpy/>`_ python package in order to use this calculator.  The
simplest way to install the former is to use your operating system's native package
manager to install the 'openkim-models' package, which will install both the KIM API and
a snapshot of binaries of all of the current models housed in the OpenKIM repository (see
`<https://openkim.org/doc/usage/obtaining-models>`_ for instructions).  Otherwise, the
'kim-api' package can be installed by itself, which will not include any models beyond
the examples bundled with the KIM API.  The kimpy package can be installed from PyPI
using pip: ``pip install --user kimpy``.

As an example, suppose we want to know the potential energy predicted by the example
model "ex_model_Ar_P_Morse_07C" for an FCC argon lattice at a lattice spacing of 5.25
Angstroms.  This can be accomplished in a manner similar to how most other ASE
calculators are used, where the name of the KIM model is passed as an argument:

::

    from ase.lattice.cubic import FaceCenteredCubic
    from ase.calculators.kim.kim import KIM

    atoms = FaceCenteredCubic(symbol='Ar', latticeconstant=5.25, size=(1,1,1))
    calc = KIM("ex_model_Ar_P_Morse_07C")
    atoms.calc = calc

    energy = atoms.get_potential_energy()
    print("Potential energy: {} eV".format(energy))

To use any other KIM model you have installed, simply substitute its name as the argument
to ``KIM``.  You can browse the models available in OpenKIM for a specific element by
visiting `<https://openkim.org>`_ and clicking on it in the periodic table.  Each model
is identified by its `extended KIM ID <https://openkim.org/doc/schema/kim-ids/>`_, which
consists of a human-readable string followed by a numeric code assigned to it.  Clicking
on an individual model will display a page containing additional information about it,
including its predictions for various material properties.  Information on how to install
KIM Models can be found at `<https://openkim.org/doc/usage/obtaining-models/>`__.

See below for a more detailed explanation of this package and additional options.

.. _Implementation:

--------------
Implementation
--------------

In order to explain the structure of this package, we must first describe the two
different types of interatomic potentials in KIM: [#typesofkimcontent]_

  * **Portable Models (PMs)**
    A KIM Portable Model (PM) is an interatomic potential designed to work
    with any simulator that supports the KIM API portable model interface.

  * **Simulator Models (SMs)**
    A KIM Simulator Model (SM) is an interatomic potential designed to work
    with a single simulator.

These two types of KIM models require different calculators to work: PMs work through a
designated calculator that uses the `kimpy <https://github.com/openkim/kimpy>`__ library
(which provides a set of python bindings to the `KIM API <https://openkim.org/kim-api>`_)
in order to set up communication between ASE and the model.  This allows, for example,
the positions of the atoms and the neighbor list in ASE to be communicated to the model,
and for the energy and forces predicted by the model for that configuration to be
communicated back to ASE in a standard format.  On the other hand, SMs are a set of
simulator commands, usually accompanied by one or more parameter files, [#smobject]_ that
are passed to a calculator corresponding to the specific simulator associated with the
SM.

Because of this separation, the :mod:`ase.calculators.kim` package consists of two modules:
:git:`ase.calculators.kim.kim <ase/calculators/kim/kim.py>` and
:git:`ase.calculators.kim.kimmodel <ase/calculators/kim/kimmodel.py>`. The first of these
contains a *wrapper function* named :func:`ase.calculators.kim.kim.KIM` that takes as
input the name of a KIM model installed on your machine, automatically determines
whether it is a PM or an SM, then constructs and returns an appropriate ASE calculator.
[#getmodelsupportedspecies]_ For example, if the name of an installed PM is passed, the
``KIM`` function will (by default) initialize an instance of
:class:`ase.calculators.kim.kimmodel.KIMModelCalculator` for it and return it as its
output.  If the name of a LAMMPS_-based SM is passed, the calculator will (by default)
return an instance of the :class:`ase.calculators.lammpslib.LAMMPSlib` calculator.  The
specific calculator type returned can be controlled using the ``simulator`` argument.

.. autofunction:: ase.calculators.kim.kim.KIM

--------------
Advanced Usage
--------------

Recalling the example given in the Overview_ section at the top of this page, no
arguments are passed to the ``KIM`` function other than the name of a portable model,
ex_model_Ar_P_Morse_07C.  From the Implementation_ section, this means that the ``calc``
object returned is actually an instance of
:class:`ase.calculators.kim.kimmodel.KIMModelCalculator` and uses the neighbor list
library implemented in kimpy.  If we wanted to use ASE's internal neighbor list
mechanism, we could specify it by modifying the corresponding line to:

::

    calc = KIM("ex_model_Ar_P_Morse_07C", options={"ase_neigh": True})

If, for some reason, we want to run our portable model with the ASE LAMMPS calculator
(:class:`ase.calculators.lammpsrun.LAMMPS`), we can specify it using the ``simulator``
argument:

::

    calc = KIM("ex_model_Ar_P_Morse_07C", simulator="lammpsrun")

Using a KIM simulator model requires no additional effort.  Using the example
LAMMPS-based simulator model bundled with the KIM API,
"Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu":

::

    from ase.lattice.cubic import FaceCenteredCubic
    from ase.calculators.kim.kim import KIM

    atoms = FaceCenteredCubic(symbol='Au', latticeconstant=4.07, size=(1,1,1))
    calc = KIM("Sim_LAMMPS_LJcut_AkersonElliott_Alchemy_PbAu")
    atoms.calc = calc

    energy = atoms.get_potential_energy()
    print("Potential energy: {} eV".format(energy))

In this case, because ``simulator`` was not specified, the default behavior is that the
object ``calc`` returned is an instance of :class:`ase.calculators.lammpslib.LAMMPSlib`.

Finally, if one instantiates the calculator for a KIM Portable Model that registers its
parameters, they can be accessed or mutated using the ``get_parameters`` and
``set_parameters`` methods.  For example, to print the components of the
parameters ``epsilons`` and ``sigmas`` in the :doi:`Lennard-Jones universal
model <10.25950/962b4967>` corresponding to Mo-Mo (index 4879), Mo-S
(index 2006) and S-S (index 1980) interactions, one can do the following:

::

    model = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    calc = KIM(model)
    print(calc.get_parameters(epsilons=[4879, 2006, 1980],
                              sigmas=[4879, 2006, 1980])

This will print a dictionary whose keys are the parameter names and whose values
are lists that each contain two sublists.  The first contains the set of indices
in the parameters arrays that were requested and the second contains the
corresponding parameter values:

::

   {'epsilons': [[4879, 2006, 1980],
                 [4.47499, 4.421814057295943, 4.36927]],
    'sigmas': [[4879, 2006, 1980],
               [2.74397, 2.30743, 1.87089]]}

Or, suppose we want to set the values of the components of the ``epsilons``
parameter of the same model corresponding to Mo-Mo, Mo-S, and S-S to 5.0, 4.5,
and 4.0, respectively.  In this case, a syntax similar is used:

::

    model = "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    calc = KIM(model)
    calc.set_parameters(epsilons=[[4879, 2006, 1980],
                                  [5.0, 4.5, 4.0]])

.. note::
   For models with parameter arrays that contain many components, such as the
   Lennard-Jones universal model above, one must refer to the documentation of
   the model itself (or its model driver, if it uses one) in order to ascertain
   the meaning of the different components of the parameter arrays.


.. _LAMMPS: http://lammps.sandia.gov

.. rubric:: Footnotes

.. [#typesofkimcontent] See `here <https://openkim.org/doc/repository/kim-content/>`_ for more details about different types of content in KIM.

.. [#smobject] Simulator Models (SMs) are actually always used in the form of a binary shared object.  This shared object contains not only the commands and other metadata related to the SM within it, but also has any relevant parameter files (if any) directly compiled into it.

.. [#getmodelsupportedspecies] As a convenience, this module also contains the function :func:`ase.calculators.kim.kim.get_model_supported_species`, which retrieves a tuple of all of the species that a model can compute interactions for.
