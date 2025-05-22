.. module:: ase.calculators.mixing

Mixing Calculators
==================

The most general class when all of the calculators has their own weight and returns with the linear combination of
calculated values. It can be used when there are different calculators eg. for the different chemical environment or
during delta leaning.

.. autoclass:: ase.calculators.mixing.LinearCombinationCalculator

There are three special variants of LinearCombinationCalculator which are available for specific tasks:

.. autoclass:: ase.calculators.mixing.MixedCalculator

.. autoclass:: ase.calculators.mixing.SumCalculator

.. autoclass:: ase.calculators.mixing.AverageCalculator

