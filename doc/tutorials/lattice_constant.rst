.. _lattice_constant:

=========================================================
Finding lattice constants using EOS and the stress tensor
=========================================================

.. seealso::

   :ref:`eos`.


HCP
===

Let's try to find the `a` and `c` lattice constants for HCP nickel
using the :mod:`EMT <ase.calculators.emt>` potential.

First, we make a good initial guess for `a` and `c` using the FCC nearest
neighbor distance and the ideal `c/a` ratio:

.. literalinclude:: lattice_constant.py
   :start-after: creates:
   :end-before: literalinclude division line 1

and create a trajectory for the results:

.. literalinclude:: lattice_constant.py
   :start-after: literalinclude division line 1
   :end-before: literalinclude division line 2

Finally, we do the 9 calculations (three values for `a` and three for `c`):

.. literalinclude:: lattice_constant.py
   :start-after: literalinclude division line 2
   :end-before: literalinclude division line 3


Analysis
--------

Now, we need to extract the data from the trajectory.  Try this:

>>> from ase.build import bulk
>>> ni = bulk('Ni', 'hcp', a=2.5, c=4.0)
>>> ni.cell
array([[ 2.5  ,  0.   ,  0.   ],
       [-1.25 ,  2.165,  0.   ],
       [ 0.   ,  0.   ,  4.   ]])

So, we can get `a` and `c` from ``ni.cell[0, 0]`` and ``ni.cell[2,
2]``:

.. literalinclude:: lattice_constant.py
   :start-after: literalinclude division line 3
   :end-before: literalinclude division line 4

We fit the energy to this expression:

.. math:: p_0 + p_1 a + p_2 c + p_3 a^2 + p_4 ac + p_5 c^2

The best fit is found like this:

.. literalinclude:: lattice_constant.py
   :start-after: literalinclude division line 4
   :end-before: literalinclude division line 5

and we can find the minimum like this:

.. literalinclude:: lattice_constant.py
   :start-after: literalinclude division line 5

Results:

.. csv-table::
   :file: lattice_constant.csv
   :header: a, c


Using the stress tensor
=======================

One can also use the stress tensor to optimize the unit cell. For this we cannot use the EMT calculator.::

    from ase.optimize import BFGS
    from ase.constraints import StrainFilter
    from gpaw import GPAW, PW
    ni = bulk('Ni', 'hcp', a=a0,c=c0)
    calc = GPAW(mode=PW(200),xc='LDA',txt='Ni.out')
    ni.calc = calc
    sf = StrainFilter(ni)
    opt = BFGS(sf)
    opt.run(0.005)

If you want the optimization path in a trajectory, add these lines
before calling the ``run()`` method::

    traj = Trajectory('path.traj', 'w', ni)
    opt.attach(traj)
