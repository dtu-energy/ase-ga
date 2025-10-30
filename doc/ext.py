from ase.utils.sphinx import mol_role
from sphinx.application import Sphinx


def setup(app: Sphinx):
    app.add_role('mol', mol_role)
