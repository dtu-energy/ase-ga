.. A new scriv changelog fragment.
..
.. Uncomment the section that is right (remove the leading dots).
.. For top level release notes, leave all the headers commented out.
..
.. I/O
.. ---
..
.. - A bullet item for the I/O category.
..
.. Calculators
.. -----------
..
.. - A bullet item for the Calculators category.
..
.. Optimizers
.. ----------
..
.. - A bullet item for the Optimizers category.
..
.. Molecular dynamics
.. ------------------
..
.. - A bullet item for the Molecular dynamics category.
..
.. GUI
.. ---
..
.. - A bullet item for the GUI category.
..
.. Development
.. -----------
..
.. - A bullet item for the Development category.
..
.. Other changes
.. -------------
..
.. - A bullet item for the Other changes category.
..
Bugfixes
--------

- Enable :func:`ase.io.magres.read_magres` to handle cases from CASTEP < 23 where indices and labels are "munged" together if the index exceeds 99. If an index exceeds 999 the situation remains ambiguous and an error will be raised. (:mr:`3530`)
