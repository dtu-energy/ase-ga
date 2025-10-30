import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', 'src').resolve()))
sys.path.append('.')

extensions = [
    'ext',  # mol role from ase
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',  # one page per object
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
]

extlinks = {
    'doi': ('https://doi.org/%s', 'doi: %s'),
    'arxiv': ('https://arxiv.org/abs/%s', 'arXiv: %s'),
}

source_suffix = '.rst'
master_doc = 'index'
project = 'ASE_GA'
# author = 'ASE developers'
# copyright = f'{datetime.date.today().year}, ASE-developers'
# exclude_patterns = ['build']
# default_role = 'math'
pygments_style = 'sphinx'
autoclass_content = 'both'
autosummary_generate = True

html_theme = 'sphinx_book_theme'
# html_logo = 'static/ase256.png'
# html_favicon = 'static/ase.ico'
# html_static_path = ['static']
html_last_updated_fmt = '%a, %d %b %Y %H:%M:%S'

html_theme_options = {
    'github_url': 'https://github.com/dtu-energy/ase-ga',
    'primary_sidebar_end': ['indices.html'],
}

# Don't sort members alphabetically
autodoc_member_order = 'bysource'

# Don't evaluate default arguments
autodoc_preserve_defaults = True

intersphinx_mapping = {
    'ase': ('https://ase-lib.org/', None),
    'python': ('https://docs.python.org/3', None),
}

# # Avoid GUI windows during doctest:
# doctest_global_setup = """
# import numpy as np
# import ase.visualize as visualize
# from ase import Atoms
# visualize.view = lambda atoms: None
# Atoms.edit = lambda self: None
# """
