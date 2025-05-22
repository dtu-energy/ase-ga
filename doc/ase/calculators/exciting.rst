.. module:: ase.calculators.exciting

========
exciting
========

.. image:: ../../static/exciting.png

Introduction
============

``exciting`` is a full-potential *all-electron*
density-functional-theory (DFT) package based on the
linearized augmented planewave (LAPW) method. It can be
applied to all kinds of materials, irrespective of the atomic species
involved, and also allows for the investigation of the core
region. The website is http://exciting-code.org/

``exciting``'s implementation in ASE requires excitingtools which is a PyPI
package that helps with the writing and reading of input/output files.
https://pypi.org/project/excitingtools/

Currently, use of ``exciting`` is limited to ground state properties.

The ExcitingGroundStateCalculator is initialized with:

* An ExcitingRunner object where the executable path is specified. If you're only planning on writing input
  files and reading output files you can feed an empty string as the executable
  path to the ExcitingRunner option.
* Ground state input options as a dictionary such as ngridk (k-points grid),
  rgkmax which dictates the muffin tin to planewave cutoff. More attributes can
  be found here: http://exciting-code.org/ref:groundstate
* The directory where to run the calculation.
* Species path directory where information about settings to use in terms of
  for example the muffin tin defaults for each species can be found.
* Optional title that get's given to the exciting ground state input xml file.

The calculator translates the ground state input options dictionary given into
the exciting input XML file format, aptly named input.xml. Note, for running a
simulation the $EXCITINGROOT environmental variable should be set: details at
http://exciting-code.org/tutorials-boron

.. literalinclude:: exciting.py

Here's an example of a ground state input xml file created using the
ExcitingGroundStateCalculator. Note, if you're not interested in running A
simulation and simply just writing the input file you can use the
ExcitingGroundStateTemplate class' write_input() method.

.. highlight:: xml

::

    <?xml version='1.0' encoding='UTF-8'?>
    <input xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://xml.exciting-code.org/excitinginput.xsd">
      <title>N3O</title>
      <structure speciespath="$EXCITINGROOT/species/" autormt="false">
        <crystal>
          <basevect>1.88972595820018 0.00000000000000 0.00000000000000</basevect>
          <basevect>0.00000000000000 1.88972595820018 0.00000000000000</basevect>
          <basevect>0.00000000000000 0.00000000000000 1.88972595820018</basevect>
        </crystal>
        <species chemicalSymbol="N" speciesfile="N.xml">
          <atom coord="0.00000000000000 0.00000000000000 0.00000000000000"/>
          <atom coord="0.00000000000000 0.00000000000000 0.00000000000000"/>
          <atom coord="0.00000000000000 0.00000000000000 0.00000000000000"/>
        </species>
        <species chemicalSymbol="O" speciesfile="O.xml">
          <atom coord="0.50000000000000 0.50000000000000 0.50000000000000"/>
        </species>
      </structure>
      <relax/>
      <groundstate tforce="true" ngridk="1 2 3"/>
      <properties>
        <dos/>
        <bandstructure>
          <plot1d>
            <path steps="100">
              <point coord="0.75000   0.50000   0.25000" label="W"/>
              <point coord="0.50000   0.50000   0.50000" label="L"/>
              <point coord="0.00000   0.00000   0.00000" label="G"/>
              <point coord="0.50000   0.50000   0.00000" label="X"/>
              <point coord="0.75000   0.50000   0.25000" label="W"/>
              <point coord="0.75000   0.37500   0.37500" label="K"/>
            </path>
          </plot1d>
        </bandstructure>
      </properties>
    </input>

The translation follows the following rules:
String values are translated to attributes. Nested dictionaries are translated
to sub elements. A list of dictionaries is translated to a list of sub elements
named after the key of which the list is the value. The special key "text()"
results in text content of the enclosing tag.


Muffin Tin Radius
=================

Sometimes it is necessary to specify a fixed muffin tin radius different from
the default. The muffin tin radii can be set by adding a custom array to the
atoms object with the name "rmt":


.. highlight:: python

::

    atoms.new_array('rmt', np.array([-1.0, -1.0, 2.3, 2.0] * Bohr))


Each entry corresponds to one atom. If the rmt value is negative, the default
value is used. This array is correctly updated if the atoms are added or removed.

Exciting Calculator Class
=========================

.. autoclass:: ase.calculators.exciting.exciting.ExcitingGroundStateCalculator




