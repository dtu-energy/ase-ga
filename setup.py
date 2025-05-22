#!/usr/bin/env python3

# Copyright (C) 2007-2017  CAMd
# Please see the accompanying LICENSE file for further information.

import os
from glob import glob
from os.path import join

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


class build_py(_build_py):
    """Custom command to build translations."""

    def __init__(self, *args, **kwargs):
        _build_py.__init__(self, *args, **kwargs)
        # Keep list of files to appease bdist_rpm.  We have to keep track of
        # all the installed files for no particular reason.
        self.mofiles = []

    def run(self):
        """Compile translation files (requires gettext)."""
        _build_py.run(self)
        msgfmt = 'msgfmt'
        status = os.system(msgfmt + ' -V')
        if status == 0:
            for pofile in sorted(glob('ase/gui/po/*/LC_MESSAGES/ag.po')):
                dirname = join(self.build_lib, os.path.dirname(pofile))
                if not os.path.isdir(dirname):
                    os.makedirs(dirname)
                mofile = join(dirname, 'ag.mo')
                print()
                print(f'Compile {pofile}')
                status = os.system(
                    '%s -cv %s --output-file=%s 2>&1' % (msgfmt, pofile, mofile)
                )
                assert status == 0, 'msgfmt failed!'
                self.mofiles.append(mofile)

    def get_outputs(self, *args, **kwargs):
        return _build_py.get_outputs(self, *args, **kwargs) + self.mofiles


setup(cmdclass={'build_py': build_py})
