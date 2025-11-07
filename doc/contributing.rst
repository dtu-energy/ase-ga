==============
 Contributing
==============

Development takes place on github https://github.com/dtu-energy/ase-ga/

Release
=======

- When ready for a release bump the version to next stable version ``uv version --bump stable``

- Push that change as the final commit of that release

- Create a release on github and give it a tag called ``v<VERSION>`` e.g. ``v1.2.3``

  A github action takes care of publishing that version to PyPI.

- Then bump the version again ``uv version --bump beta`` and push that commit, signifying that new changes are ready to be added.
