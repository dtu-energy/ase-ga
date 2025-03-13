.. _writing_documentation_ase:

=====================
Writing documentation
=====================

We use the Sphinx_ tool to generate the documentation.  The documentation is
stored on GitLab as text files in the :git:`doc` directory using the
reStructuredText_ markup language.

.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _Sphinx: http://www.sphinx-doc.org/en/master/


Installing Docutils and Sphinx
==============================

If you do:

.. code-block:: console

    $ pip install sphinx_rtd_theme --user

and add ``~/.local/bin`` to you :envvar:`PATH` environment variable, then
you should be ready to go.  You may need the following installed, but they
are not required: scipy, matplotlib, povray, dvipng, pdflatex, bibtex,
AUCTex, fontconfig, convert (ImageMagick).


.. _using_sphinx:

Using Sphinx
============

First, you should take a look at the documentation for Sphinx_ and
reStructuredText_.

If you don't already have your own copy of the ASE package, then read
:ref:`here <contribute>` how to get everything set up.

Then :command:`cd` to the :file:`doc` directory and build the HTML pages:

.. code-block:: console

  $ cd ~/ase/doc
  $ make

This might take a long time the first time you do it.

.. Note::

   Make sure that you build the Sphinx documentation using the
   corresponding ASE version by setting the environment variables
   :envvar:`PYTHONPATH` and :envvar:`PATH`.

Create a branch for your work, make your changes to the ``.rst`` files, run
:command:`make` again, check the results and if things
look ok, create a *merge request*:

.. code-block:: console

    $ git checkout -b fixdoc
    $ idle index.rst
    $ make
    $ make browse
    $ git commit -am "..."
    $ git push -u origin fixdoc


Extensions to Sphinx
====================

.. highlight:: rest

We have a couple of extensions to Sphinx:

**:mol:**

   Use ``:mol:`CH_3OH``` to get :mol:`CH_3OH`.

**:git:**

   A role for creating a link to a file on GitLab.  If you write
   ``:git:`ase/atoms.py```, you
   will get: :git:`ase/atoms.py`.

**:math:**

   This role is for inline LaTeX-style math.  Example:
   ``:math:`\sin(x_n^2)``` gives you :math:`\sin(x_n^2)`.  This role
   is actually the default for ASE's documentation, so you should leave
   out the ``:math:`` part like here: ```\sin(x_n^2)```.

**.. math::**

   Write displayed LaTeX-style math.  Example::

     .. math:: \frac{1}{1+x^2}

   gives you:

   .. math:: \frac{1}{1+x^2}


.. _generated:

Running Python code to create figures
=====================================

If you want to include a picture in your page, *you should not* check
in the png-file to our Git repositoy!  Instead, you should check in
the Python script you used to generate the picture (you can also
generate csv-files or pdf-files like this).  The first line of the
script should look like this::

    # creates: fig1.png, fig2.png, table1.csv

Sphinx will run the script and generate the files that you can
then use in your rst-file.  Examples:

* :ref:`eos`.  Source: :git:`doc/tutorials/eos/eos.py`,
  :git:`doc/tutorials/eos/eos.rst`
* :ref:`lattice_constant`.  Source: :git:`doc/tutorials/lattice_constant.py`,
  :git:`doc/tutorials/lattice_constant.rst`


reStructedText in emacs
=======================

.. highlight:: common-lisp

For people using emacs, the `reStructuredText extension`_ is highly
recommended. The intallation procedure is described in the top of the
file, but for most people, it is enough to place it in your emacs
load-path (typically ``.emacs.d/``) and add the lines::

  (add-to-list 'load-path "~/.emacs.d")
  (require 'rst)

somewhere in your ``.emacs`` file.

To make the mode auto load for relevant file extension, you can write
something like::

  (setq auto-mode-alist
        (append '(("\\.rst$" . rst-mode)
                  ("\\.rest$" . rst-mode)) auto-mode-alist))

In your ``.emacs`` file.

.. _reStructuredText extension: http://docutils.sourceforge.net/
                                tools/editors/emacs/rst.el
