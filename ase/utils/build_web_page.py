"""Build ASE's web-page."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

cmds = """\
. /etc/bashrc
module load Tkinter
python3 -m venv venv
. venv/bin/activate
pip install -U pip
pip install "sphinx<6.0"  # search broken in sphinx-6
pip install sphinx-rtd-theme pillow
git clone http://gitlab.com/ase/ase.git
cd ase
pip install .
python setup.py sdist
cd doc
make
mv build/html ase-web-page"""


def build():
    root = Path('/tmp/ase-docs')
    if root.is_dir():
        sys.exit('Locked')
    root.mkdir()
    os.chdir(root)
    cmds2 = ' && '.join(line.split('#')[0] for line in cmds.splitlines())
    p = subprocess.run(cmds2, shell=True, check=False)
    if p.returncode == 0:
        status = 'ok'
    else:
        print('FAILED!', file=sys.stdout)
        status = 'error'
    f = root.with_name(f'ase-docs-{status}')
    if f.is_dir():
        shutil.rmtree(f)
    root.rename(f)
    return status


def build_all():
    assert build() == 'ok'
    tar = next(
        Path('/tmp/ase-docs-ok/ase/dist/').glob('ase-*.tar.gz'))
    webpage = Path('/tmp/ase-docs-ok/ase/doc/ase-web-page')
    home = Path.home() / 'web-pages'
    cmds = ' && '.join(
        [f'cp {tar} {webpage}',
         f'find {webpage} -name install.html | '
         f'xargs sed -i s/snapshot.tar.gz/{tar.name}/g',
         f'cd {webpage.parent}',
         'tar -czf ase-web-page.tar.gz ase-web-page',
         f'cp ase-web-page.tar.gz {home}'])
    subprocess.run(cmds, shell=True, check=True)


if __name__ == '__main__':
    build_all()
