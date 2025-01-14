.. _coding conventions:

==================
Coding Conventions
==================

The code must be compatible with the oldest supported version of python
as given on the :ref:`download_and_install` page.

Please run :ref:`ruff check <stylecheck>` on your code.

The rules are almost identical
to those used by the `Docutils project`_:

Contributed code will not be refused merely because it does not
strictly adhere to these conditions; as long as it's internally
consistent, clean, and correct, it probably will be accepted.  But
don't be surprised if the "offending" code gets fiddled over time to
conform to these conventions.

The project follows the generic coding conventions as
specified in the `Style Guide for Python Code`_ and `Docstring
Conventions`_ PEPs, clarified and extended as follows:

* Do not use "``*``" imports such as ``from module import *``.  Instead,
  list imports explicitly.

* Use 4 spaces per indentation level.  No tabs.

* Read the *Whitespace in Expressions and Statements*
  section of PEP8_.

* Avoid `trailing whitespaces`_.

* No one-liner compound statements (i.e., no ``if x: return``: use two
  lines).

* Maximum line length is 78 characters.

* Use "StudlyCaps" for class names.

* Use "lowercase" or "lowercase_with_underscores" for function,
  method, and variable names.  For short names,
  joined lowercase may be used (e.g. "tagname").  Choose what is most
  readable.

* No single-character variable names, except indices in loops
  that encompass a very small number of lines
  (``for i in range(5): ...``).

* Avoid lambda expressions.  Use named functions instead.

* Avoid functional constructs (filter, map, etc.).  Use list
  comprehensions instead.

* Use ``'single quotes'`` for string literals, and ``"""triple double
  quotes"""`` for :term:`docstring`\ s.  Double quotes are OK for
  something like ``"don't"``.

.. _Style Guide for Python Code:
.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _Docstring Conventions: https://www.python.org/dev/peps/pep-0257/
.. _Docutils project: http://docutils.sourceforge.net/docs/dev/policies.html
                      #python-coding-conventions
.. _trailing whitespaces: http://www.gnu.org/software/emacs/manual/html_node/
                          emacs/Useless-Whitespace.html

.. attention::

   Thus spake the Lord: Thou shalt indent with four spaces. No more, no less.
   Four shall be the number of spaces thou shalt indent, and the number of thy
   indenting shall be four. Eight shalt thou not indent, nor either indent thou
   two, excepting that thou then proceed to four. Tabs are right out.

                                          Georg Brandl


General advice
==============

 * Get rid of as many ``break`` and ``continue`` statements as possible.

 * Write short functions.  All functions should fit within a standard screen.

 * Use descriptive variable names.

Writing documentation in the code
=================================

ASE follows the NumPy/SciPy convention for docstrings:

  https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard


.. _stylecheck:

Run ruff on your code
=======================

It's a good idea to run `ruff <https://docs.astral.sh/ruff/>`_
on your code (or use a text editor that does it automatically)::

    $ ruff check --fix filename.py
    $ ruff format filename.py

.. _autopep8py:

Run autopep8.py on your code
============================

Another method of enforcing PEP8_ is using a tool such as
`autopep8.py <https://github.com/hhatto/autopep8>`_. These tools tend to be
very effective at cleaning up code, but should be used carefully and code
should be retested after cleaning it. Try::

  $ autopep8.py --help

.. attention::

   There is a common issue with pep8 where spaces are added around the power
   operator.  Code such as "x**2" should not be changed to "x ** 2".  This
   issue is not fixed in pep8 as of the time of this writing, but a small
   `change <http://listserv.fysik.dtu.dk/pipermail/gpaw-developers/
   2014-October/005075.html>`_ to autopep8 has been effective to prevent
   this change.
