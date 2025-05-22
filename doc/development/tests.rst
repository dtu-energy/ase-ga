.. module:: ase.test

================
Testing the code
================

All additions and modifications to ASE should be tested.

.. index:: testase

Test scripts should be put in the :git:`ase/test` directory.
Run all tests with::

  ase test

This requires installing ``pytest``, ``pytest-mock`` and ``pytest-xdist``.
You can install these dependencies automatically by running::

  pip install ase[test]

Tests for some calculators require data files in
https://gitlab.com/ase/ase-datafiles.
You can enable these tests by running::

  pip install --user --upgrade git+https://gitlab.com/ase/ase-datafiles.git

See ``ase test --help`` for more information.

You can also run ``pytest`` directly from within the ``ase.test`` directory.

.. important::

  When you fix a bug, add a test to the test suite checking that it is
  truly fixed.  Bugs sometimes come back, do not give it a second
  chance!


How to add a test
=================

Create a module somewhere under ``ase.test``.  Make sure its name
starts with ``test_``.  Inside the module, each test should be a
function whose name starts with ``test_``.  This ensures that pytest
finds the test.  Use ``ase test --list`` to see which tests it will
find.

You may note that many tests do not follow these rules.
These are older tests.  We expect to port them one day.

How to fail successfully
========================

The test suite provided by :mod:`ase.test` automatically runs all test
scripts in the :git:`ase/test` directory and summarizes the results.

If a test script causes an exception to be thrown, or otherwise terminates
in an unexpected way, it will show up in this summary. This is the most
effective way of raising awareness about emerging conflicts and bugs during
the development cycle of the latest revision.


Remember, great tests should serve a dual purpose:

**Working interface**
    To ensure that the :term:`class`'es and :term:`method`'s in ASE are
    functional and provide the expected interface. Empirically speaking, code
    which is not covered by a test script tends to stop working over time.

**Replicable results**
    Even if a calculation makes it to the end without crashing, you can never
    be too sure that the numerical results are consistent. Don't just assume
    they are, :func:`assert` it!

.. function:: assert(expression)

    Raises an ``AssertionError`` if the ``expression`` does not
    evaluate to ``True``.



Example::

  from ase import molecule

  def test_c60():
      atoms = molecule('C60')
      atoms.center(vacuum=4.0)
      result = atoms.get_positions().mean(axis=0)
      expected = 0.5*atoms.get_cell().diagonal()
      tolerance = 1e-4
      assert (abs(result - expected) < tolerance).all()


To run the same test with different inputs, use pytest fixtures.
For example::

  @pytest.mark.parametrize('parameter', [0.1, 0.3, 0.7])
  def test_something(parameter):
      # setup atoms here...
      atoms.set_something(parameter)
      # calculations here...
      assert everything_is_going_to_be_alright
