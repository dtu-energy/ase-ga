Resonant and non-resonant Raman spectra
=======================================

.. note::
   Raman spectra can also be obtained via :mod:`~ase.calculators.siesta` calculator.

Raman spectra can be calculated in various approximations [1]_.
While the examples below are using GPAW_ explicitly,
the modules are intended to work with other calculators also.
The strategy is to calculate vibrational properties first and
obtain the spectra from these later.

---------------------------------
1. Finite difference calculations
---------------------------------

1a. Forces and excitations
--------------------------

The basis for all spectra are finite difference calculations
for forces and excited states of our system.
These can be performed using the
:class:`~ase.vibrations.resonant_raman.ResonantRamanCalculator`.

.. autoclass:: ase.vibrations.resonant_raman.ResonantRamanCalculator


1b. More accurate forces
------------------------

It is possible to do a vibrational also with a more accurate calculator
(or more accurate settings for the forces) using
the :class:`~ase.vibrations.Vibrations` or  :class:`~ase.vibrations.Infrared`
modules.

In the example of molecular hydrogen with GPAW_ this is

.. literalinclude:: H2_ir.py

This produces a calculation with rather accurate forces in order
to get the Hessian and thus the vibrational frequencies as well
as Eigenstates correctly.

In the next step we perform a finite difference optical calculation
with less accuracy,
where the optical spectra are evaluated using TDDFT

.. literalinclude:: H2_optical.py
		    

1c. Overlaps
------------

Albrecht B+C terms need wave function overlaps between equilibrium and
displaced structures. These are assumed to be
calculated in the form

.. math::

  o_{ij} = \int d\vec{r} \; \phi_i^{{\rm disp},*}(\vec{r})
  \phi_j^{{\rm eq}}(\vec{r})
   
where `\phi_j^{{\rm eq}}` is an orbital at equilibrium position
and `\phi_i^{\rm disp}` is an orbital at displaced position.

The ``H2MorseExcitedStatesCalculator`` has a function ``overlap()`` for this.
We therfore write data including the overlap as

.. literalinclude:: H2Morse_calc_overlap.py

In GPAW this is implemented in ``Overlap``
(approximated by pseudo-wavefunction overlaps) and can be triggered
in ``ResonantRamanCalculator`` by


.. literalinclude:: H2_optical_overlap.py


2. Analysis of the results
--------------------------

We assume that the steps above were performed and are able to analyse the
results in different approximations.

Placzek
```````

The most popular form is the Placzeck approximation that is present in
two implementations. The simplest is the direct evaluation from
derivatives of the frequency dependent polarizability::

  from ase.calculators.h2morse import (H2Morse,
                                       H2MorseExcitedStates)
  from ase.vibrations.placzek import Placzek

  photonenergy = 7.5  # eV
  pz = Placzek(H2Morse(), H2MorseExcitedStates)
  pz.summary(photonenergy)


The second implementation evaluates the derivatives differently, allowing
for more analysis::

  import pylab as plt
  from ase.calculators.h2morse import (H2Morse,
                                       H2MorseExcitedStates)
  from ase.vibrations.placzek import Profeta
  
  photonenergy = 7.5  # eV
  pr = Profeta(H2Morse(), H2MorseExcitedStates, approximation='Placzek')
  x, y = pr.get_spectrum(photonenergy, start=4000, end=5000, type='Lorentzian')
  plt.plot(x, y)
  plt.show()

Both implementations should lead to the same spectrum.

``Profeta`` splits the spectra in two contributions that can be accessed as
``approximation='P-P'`` and ``approximation='Profeta'``, respectively.
Their sum should give ``approximation='Placzek'``.
See more details in [1]_.

Albrecht
````````

The more accurate Albrecht approximations partly need overlaps
to be present. We therefore have to invoke the ``Albrecht`` object as::

  from ase.calculators.h2morse import (H2Morse,
                                       H2MorseExcitedStates)
  from ase.vibrations.albrecht import Albrecht
  
  photonenergy = 7.5  # eV
  al = Albrecht(H2Morse(), H2MorseExcitedStates, approximation='Albrecht', overlap=True)
  x, y = al.get_spectrum(photonenergy, start=4000, end=5000, type='Lorentzian')

``Albrecht`` splits the spectra in two contributions that can be accessed as
``approximation='Albrecht A'`` and ``approximation='Albrecht BC'``,
respectively.
Their sum should give ``approximation='Albrecht'``.
See more details in [1]_.

  
.. _GPAW: https://gpaw.readthedocs.io/
  
.. [1] :doi:`Ab-initio wave-length dependent Raman spectra: Placzek approximation and beyond <10.1021/acs.jctc.9b00584>`

