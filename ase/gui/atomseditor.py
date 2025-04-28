import ase.gui.ui as ui
from ase.gui.i18n import _


class AtomsEditor:
    def __init__(self, gui):
        gui.register_vulnerable(self)

        win = ui.Window(_('Edit atoms'), wmtype='utility')

        tree = ui.ttk.Treeview(win.win, selectmode='browse')
        tree.pack(side='left', fill='y')
        bar = ui.ttk.Scrollbar(win.win, orient='vertical', command=tree.yview)
        tree.configure(yscrollcommand=bar.set)
        bar.pack(side='right', fill='y')

        columns = ['symbol', *'xyz']
        tree['columns'] = columns

        tree.column('#0', width=40)
        tree.heading('#0', text=_('id'))

        tree.column('symbol', width=60)
        tree.heading('symbol', text=_('symbol'))

        for name in 'xyz':
            tree.column(name, width=80)
            tree.heading(name, text=name)

        atoms = gui.atoms
        for i in range(len(atoms)):
            values = [atoms.symbols[i], *atoms.positions[i]]
            tree.insert('', 'end', text=i, values=values)
