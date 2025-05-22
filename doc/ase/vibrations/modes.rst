.. module:: ase.vibrations

Vibrational modes
=================

You can calculate the vibrational modes of an
:class:`~ase.Atoms` object in the harmonic approximation using
the :class:`Vibrations`.

.. autoclass:: Vibrations
   :members:

Example
-------

The example of a water molecule in the EAM potential

.. literalinclude:: H2O_EMT.py

where the output

.. include:: H2O_EMT_summary.txt
   :literal:

shows 3 meaningful vibrations (the last 3
with highest energies.)

These vibrations can be viewed in ``ase gui``
either by writing them out as a "movie"::

  vib.write_mode(-1)

which writes out the file ``vib.8.traj``
The vibrations can also be encoded as forces::

  vib.show_as_force(8)

which opens ``ase gui`` automatically and the forces point
into directions of the movement of the atoms. 


Old calculations
----------------

The output format of vibrational calculations was changed from ``pickle``
to ``json``. There is a tool to convert old ``pickle``-files::

  > python3 -m ase.vibrations.pickle2json mydirectory/vib.*.pckl

Vibrational Data
----------------
Vibrational data is stored inside the :class:`~ase.vibrations.data.VibrationsData` class.

.. autoclass:: VibrationsData
   :members:
