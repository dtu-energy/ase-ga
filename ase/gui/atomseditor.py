from dataclasses import dataclass
from functools import partial
from typing import Callable

import ase.gui.ui as ui
from ase.gui.i18n import _


@dataclass
class Column:
    name: str
    displayname: str
    widget_width: int
    getvalue: Callable
    setvalue: Callable


class AtomsEditor:
    def __init__(self, gui):
        gui.register_vulnerable(self)

        win = ui.Window(_('Edit atoms'), wmtype='utility')

        tree = ui.ttk.Treeview(win.win, selectmode='extended')
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

        tree.bind('<Return>', self.pressed_return)
        tree.bind('<Double-1>', self.doubleclick)

        self.tree = tree

    def pressed_return(self, event):
        print(event)

    def doubleclick(self, event):
        print(event)

        # what row and column was clicked on
        row_id = self.tree.identify_row(event.y)
        column_id = self.tree.identify_column(event.x)
        print('click', row_id, column_id)

        assert column_id[0] == '#'
        column = int(column_id[1:])


        x, y, width, height = self.tree.bbox(row_id, column_id)

        item = self.tree.item(row_id, 'values')
        print(item)
        if column == 0:
            return

        print('column', column)
        content = item[column - 1]
        print(content)

        entry = ui.ttk.Entry(self.tree)
        entry.insert(0, content)
        entry.focus_force()
        entry.selection_range(0, 'end')

        def apply_change(event):
            print('apply', event, entry, entry.get())
            entry.destroy()

        entry.bind('<Return>', apply_change)
        entry.bind('<Escape>', lambda *args: entry.destroy())

        # entry.place(row_id, column - 1, 'hello')
        entry.place(x=x, y=y, width=width, height=height)

