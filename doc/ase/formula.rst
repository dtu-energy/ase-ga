.. module:: ase.formula

Chemical formula type
=====================

The :class:`~Formula` type collects all formula manipulation functionality
in one place (see examples below).

The following string formats are supported:

.. csv-table::
   :file: formats.csv
   :header-rows: 1

.. autoclass:: Formula
    :members:
    :member-order: bysource
    :special-members: __format__, __contains__, __getitem__,
       __eq__, __divmod__, __add__, __mul__, __len__
