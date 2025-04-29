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

        tree.column('#0', width=40)
        tree.heading('#0', text=_('id'))

        def get_symbol(atoms, i):
            return atoms.symbols[i]

        def set_symbol(atoms, i, value):
            from ase.data import atomic_numbers

            if value not in atomic_numbers:
                pass  # Display error?
            atoms.symbols[i] = value

        self.gui = gui
        self.tree = tree
        self.columns = []

        self.add_column(
            'symbol',
            _('symbol'),
            60,
            lambda atoms, i: atoms.symbols[i],
            set_symbol,
        )

        self.update_table_from_atoms()

        tree.bind('<Return>', self.pressed_return)
        tree.bind('<Double-1>', self.doubleclick)

    @property
    def atoms(self):
        return self.gui.atoms

    def update_table_from_atoms(self):
        for i in range(len(self.atoms)):
            values = self.get_row_values(i)
            self.tree.insert('', 'end', text=i, values=values)

    def get_row_values(self, i):
        return [column.getvalue(self.gui.atoms, i) for column in self.columns]

    def add_column(self, name, heading, width, getvalue, setvalue):
        column = Column(name, heading, width, getvalue, setvalue)

        names = [column.name for column in self.columns]
        names.append(name)

        self.tree['columns'] = names
        self.tree.column(name, width=width)
        self.tree.heading(name, text=heading)
        self.columns.append(column)

    def pressed_return(self, event):
        print(event)

    def doubleclick(self, event):
        print(event)

        row_id = self.tree.identify_row(event.y)
        column_id = self.tree.identify_column(event.x)

        assert column_id[0] == '#'
        assert row_id.startswith('I')

        column_no = int(column_id[1:]) - 1
        if column_no == -1:
            return  # This is the ID column.
        row_no = int(row_id[1:].rstrip('0')) - 1

        assert 0 <= column_no < len(self.columns)
        assert 0 <= row_no < len(self.atoms)

        x, y, width, height = self.tree.bbox(row_id, column_id)

        item = self.tree.item(row_id, 'values')

        content = item[column_no]

        entry = ui.ttk.Entry(self.tree)
        entry.insert(0, content)
        entry.focus_force()
        entry.selection_range(0, 'end')

        def apply_change(event):
            print('apply', event, entry, entry.get())
            column = self.columns[column_no]
            print('Column', column)
            value = entry.get()
            column.setvalue(self.atoms, row_no, value)
            self.tree.set(row_id, column_id, value=value)
            self.gui.set_frame()
            entry.destroy()

        entry.bind('<Return>', apply_change)
        entry.bind('<Escape>', lambda *args: entry.destroy())

        entry.place(x=x, y=y, width=width, height=height)
