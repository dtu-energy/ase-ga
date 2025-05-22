==============
Bader Analysis
==============

An example of Bader analysis: `The water molecule`_.

.. _The water molecule: https://gpaw.readthedocs.io/tutorialsexercises/
                        wavefunctions/charge/charge.html#bader

You can attach the output charges from the bader program to the atoms
for further processing::

    from ase.io.bader import attach_charges

    # the next two lines are equivalent (only one needed)
    attach_charges(atoms)
    attach_charges(atoms, 'ACF.dat')

    for atom in atoms:
        print('Atom', atom.symbol, 'Bader charge', atom.charge)
