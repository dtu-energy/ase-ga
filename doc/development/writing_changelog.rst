.. _writing_changelog:

=================
Writing changelog
=================

Since ASE 3.24.0 (:mr:`3572`),
we recommend using |scriv|_ to add your changes in the changelog rather than
updating ``doc/releasenotes.rst`` directly to avoid merge conflicts on it.

.. |scriv| replace:: ``scriv``
.. _scriv: https://scriv.readthedocs.io/

Installing ``scriv``
====================

If you have not installed ``scriv``, you should first install it, e.g., as:

.. code-block:: console

    $ pip install scriv

Using ``scriv``
===============

Once you have made changes on ASE, you run ``scriv``, e.g., as:

.. code-block:: console

    $ scriv create --add

It makes a file like ``20250101_000000_john_doe_my_change.rst``
in the ``changelog.d`` directory.
You can also rename the file under the rule ``<timestamp>_<subject>.rst``,
where ``<timestamp>`` should be at least ``YYYYMMDD``.

You then uncomment the relevant section, add some notes about the change,
and commit the updated file.

How does it work?
=================

When ASE maintainers make a new release, they will compile these files
automatically using the command ``scriv collect``.
This will put all changes given in ``changelog.d`` in ``doc/releasenotes.rst``.
