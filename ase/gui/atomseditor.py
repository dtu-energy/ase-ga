import ase.gui.ui as ui
from ase.gui.i18n import _


class AtomsEditor:
    def __init__(self, gui):
        gui.register_vulnerable(self)

        win = ui.Window(_('Edit atoms'), wmtype='utility')

        self.gui = gui
        self.win = win
