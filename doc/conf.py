import datetime

extensions = [
    'ase.utils.sphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.imgconverter',
]

extlinks = {
    'doi': ('https://doi.org/%s', 'doi: %s'),
    'arxiv': ('https://arxiv.org/abs/%s', 'arXiv: %s'),
    'mr': ('https://gitlab.com/ase/ase/-/merge_requests/%s', 'MR: !%s'),
    'issue': ('https://gitlab.com/ase/ase/-/issues/%s', 'issue: #%s'),
}
source_suffix = '.rst'
master_doc = 'index'
project = 'ASE'
copyright = f'{datetime.date.today().year}, ASE-developers'
templates_path = ['templates']
exclude_patterns = ['build']
default_role = 'math'
pygments_style = 'sphinx'
autoclass_content = 'both'
modindex_common_prefix = ['ase.']
nitpick_ignore = [
    ('envvar', 'VASP_PP_PATH'),
    ('envvar', 'ASE_ABC_COMMAND'),
    ('envvar', 'LAMMPS_COMMAND'),
    ('envvar', 'ASE_NWCHEM_COMMAND'),
    ('envvar', 'SIESTA_COMMAND'),
    ('envvar', 'SIESTA_PP_PATH'),
    ('envvar', 'VASP_SCRIPT'),
]

html_theme = 'sphinx_rtd_theme'
html_style = 'ase.css'
html_favicon = 'static/ase.ico'
html_static_path = ['static']
html_last_updated_fmt = '%a, %d %b %Y %H:%M:%S'

latex_elements = {'papersize': 'a4paper'}
latex_show_urls = 'inline'
latex_show_pagerefs = True
latex_engine = 'xelatex'
latex_documents = [
    ('index', 'ASE.tex', 'ASE', 'ASE-developers', 'howto', not True)
]

intersphinx_mapping = {
    'gpaw': ('https://gpaw.readthedocs.io', None),
    'python': ('https://docs.python.org/3.10', None),
}

# Avoid GUI windows during doctest:
doctest_global_setup = """
import numpy as np
import ase.visualize as visualize
from ase import Atoms
visualize.view = lambda atoms: None
Atoms.edit = lambda self: None
"""

autodoc_mock_imports = ['kimpy']
