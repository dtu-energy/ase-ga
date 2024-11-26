.. module:: ase.phasediagram
.. _phase diagrams:

====================================
Phase diagrams and Pourbaix diagrams
====================================

.. autoclass:: ase.phasediagram.PhaseDiagram

Here is a simple example using some made up numbers for Cu-Au alloys:

>>> from ase.phasediagram import PhaseDiagram
>>> refs = [('Cu', 0.0),
...         ('Au', 0.0),
...         ('CuAu', -0.5),
...         ('Cu2Au', -0.7),
...         ('Cu2Au', -0.2)]
>>> pd = PhaseDiagram(refs)
Species: Au, Cu
References: 5
0    Cu             0.000
1    Au             0.000
2    CuAu          -0.500
3    Cu2Au         -0.700
4    CuAu2         -0.200
Simplices: 3

The convex hull looks like this:

>>> pd.plot(show=True)

.. image:: cuau.png

.. automethod:: PhaseDiagram.plot

If you want to see what :mol:`Cu_3Au` will decompose into, you can use the
:meth:`~PhaseDiagram.decompose` method:

>>> energy, indices, coefs = pd.decompose('Cu3Au')
reference    coefficient      energy
------------------------------------
Cu                     1       0.000
Cu2Au                  1      -0.700
------------------------------------
Total energy:                 -0.700
------------------------------------
>>> print(energy, indices, coefs)
(-0.69999999999999996, array([0, 3], dtype=int32), array([ 1.,  1.]))

Alternatively, one could have used ``pd.decompose(Cu=3, Au=1)``.

.. automethod:: PhaseDiagram.decompose

Here is an example (see :download:`ktao.py`) with three components using
``plot(dims=2)`` and ``plot(dims=3)``:

.. image:: ktao-2d.png
.. image:: ktao-3d.png


Pourbaix diagrams
=================

Let's create a Pourbaix diagram for ZnO from experimental numbers.
We start by collecting the experimental formation energies of
solvated zinc-containing ions and molecules by using the :func:`solvated` function:

>>> from ase.pourbaix import Pourbaix
>>> from ase.phasediagram import solvated
>>> refs = solvated('Zn')
>>> print(refs)
[('HZnO2-(aq)', -4.801274772854441), ('ZnO2--(aq)', -4.0454382546928365), ('ZnOH+(aq)', -3.5207324675582736), ('ZnO(aq)', -2.9236086089762137), ('H2O(aq)', -2.458311658897383), ('Zn++(aq)', -1.5264168353005447), ('H+(aq)', 0.0)]

.. autofunction:: solvated

We add two solids and one more dissolved molecule to the references,
convert them to a dictionary and create a :class:`~ase.pourbaix.Pourbaix`
object using the class default arguments:

>>> refs += [('Zn', 0.0), ('ZnO', -3.323), ('ZnO2(aq)', -2.921)]
>>> pb = Pourbaix('ZnO', dict(refs))

We can determine what is the most stable phase at a potential of 1 V
and a pH of 9.0, see the corresponding chemical reaction and determine
the Pourbaix energy, i.e. the energy of the target material ZnO
relative to the most stable competing phase:

>>> energy, phase = pbx.get_pourbaix_energy(1.0, 9.0, verbose=True)
Stable phase:
ZnO + H2O  ➜  2H+ + 2e- + ZnO2(aq)
Energy: 0.560 eV

As we can see in these conditions ZnO spontaneously decomposes
into acqueous :mol:`ZnO2` and lies 560 meV above the latter species.
This chemical reaction is described by a :class:`ase.pourbaix.RedOx` object.
We can show that, in the same conditions, the reaction occurs spontaneously
releasing the opposite of the pourbaix energy, i.e. -560 meV:

>>> print(type(phase))
<class 'ase.pourbaix.RedOx'>
>>> print(phase.get_free_energy(1.0, 9.0))
-0.5595239105125918

If we repeat the evaluation at a potential of 0 V and a pH of 10,
we can see that ZnO is now the most stable phase (although barely),
and the associated Pourbaix energy is (slightly) negative.

>>> pbx.get_pourbaix_energy(0.0, 10.0, verbose=True)
ZnO is stable.
Energy: -0.033 eV

Finally, we can evaluate the complete Pourbaix diagram of ZnO
in a potential window between -2 and +2 V, and a pH window between
0 and 14:

>>> Urange = [-2, 2]
>>> pHrange = [0, 14]
>>> pbx.plot(Urange, pHrange, show=True)

.. image:: zno.png

.. autoclass:: ase.pourbaix.Pourbaix
.. autoclass:: ase.pourbaix.RedOx
