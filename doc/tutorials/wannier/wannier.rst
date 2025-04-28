.. _wannier tutorial:
    
=================================
Partly occupied Wannier Functions
=================================
In this tutorial we show how to use the ASE Wannier functions module. To run the tutorial, you will need to have GPAW installed.

1: Benzene module

First we will calculate Wannier functionsd for the benzene molecule. First we will need to run a ground state calculation to find the ground state electron density and Kohn-Sham wave functions. The result of the ground state calculation is saved to the file benzene.gpw.
.. literalinclude:: benzene.py

Now, we are ready to construct the Wannier functions.
We first make a Wannierization with 15 Wannier functions, which matches the number of occupied bands in the benzene molecule. This is specified using the nwannier parameter in the constructor. 
The call to the localize() functions attempt to localize the Wannier functions using a gradient descent optimization procedure. The resulting maximally localized Wannier functions are saved to a .cube file, which allows them to be visualized.

Next, we examine the effect of including extra degrees of freedom. These are extra Wannier functions that are included to improve the localization, but which need not describe any particular bands of the system. The extra degrees of freedom may be included by setting 'fixedstates' to a lower number than the total number of Wannier functions, given by 'nwannier'. By specifying 18 Wannier functions in total, we should therefore obtain a set of Wannier functions that are more localized than before.
.. literalinclude:: wannier_benzene.py

We can plot the projections of the Wannier functions on the energy eigenstates. The resulting plot shows how the 15 lowest bands, corresponding to the fixed states that we specified, are perfectly represented by the Wannier functions. Notice how the Wannier functions also include contributions from the unoccupied bands, but that none of these bands are perfectly represented. This is due to the optimization of the extra degrees of freedom that we wanted.
.. literalinclude:: plot_spectral_weight.py


.. literalinclude:: polyacetylene.py
.. literalinclude:: wannier_polyacetylene.py
.. literalinclude:: plot_band_structure.py
