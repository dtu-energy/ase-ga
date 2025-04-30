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
    format_value: Callable = lambda obj: str(obj)


class AtomsEditor:
    def __init__(self, gui):
        gui.register_vulnerable(self)

        win = ui.Window(_('Edit atoms'), wmtype='utility')

        tree = ui.ttk.Treeview(win.win, selectmode='extended')
        edit_entry = ui.ttk.Entry(win.win)
        edit_entry.pack(side='bottom', fill='x')#fill='x')
        tree.pack(side='left', fill='y')
        bar = ui.ttk.Scrollbar(win.win, orient='vertical', command=tree.yview)
        tree.configure(yscrollcommand=bar.set)

        tree.column('#0', width=40)
        tree.heading('#0', text=_('id'))

        bar.pack(side='right', fill='y')

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

        self.define_column(
            'symbol',
            _('symbol'),
            60,
            lambda atoms, i: atoms.symbols[i],
            set_symbol,
        )

        for c, axisname in enumerate('xyz'):
            class GetSetPos:
                def __init__(self, c):
                    self.c = c

                def set_position(self, atoms, i, value):
                    atoms.positions[i, self.c] = float(value)
                def get_position(self, atoms, i):
                    return atoms.positions[i, self.c]

            self.define_column(
                axisname,
                axisname,
                92,
                GetSetPos(c).get_position,
                GetSetPos(c).set_position,
                format_value=lambda val: f'{val:.4f}',
            )

        self.update_table_from_atoms()

        tree.bind('<Return>', self.pressed_return)
        tree.bind('<Double-1>', self.doubleclick)

        self.edit_entry = edit_entry

    # notify when atoms are changed (new_atoms())
    # be able to listen to any editing and update fields

    @property
    def atoms(self):
        return self.gui.atoms

    def update_table_from_atoms(self):
        for i in range(len(self.atoms)):
            values = self.get_row_values(i)
            self.tree.insert('', 'end', text=i, values=values)
        self.add_columns_to_widget()

    def get_row_values(self, i):
        return [column.format_value(column.getvalue(self.gui.atoms, i))
                for column in self.columns]

    def define_column(self, *args, **kwargs):
        column = Column(*args, **kwargs)
        self.columns.append(column)

    def add_columns_to_widget(self):
        self.tree['columns'] = [column.name for column in self.columns]
        for column in self.columns:
            self.tree.heading(column.name, text=column.displayname)
            self.tree.column(column.name, width=column.widget_width,
                             anchor='e')

    def pressed_return(self, event):
        print(event)

    def doubleclick(self, event):
        row_id = self.tree.identify_row(event.y)
        column_id = self.tree.identify_column(event.x)

        assert column_id.startswith('#'), repr(column_id)
        assert row_id.startswith('I'), repr(row_id)

        column_no = int(column_id[1:]) - 1
        if column_no == -1:
            return  # This is the ID column.

        # WTF why in the name of the devil does it use hexadecimal
        row_no = int(row_id[1:], base=16) - 1

        assert 0 <= column_no < len(self.columns)
        assert 0 <= row_no < len(self.atoms)

        content = self.columns[column_no].getvalue(self.atoms, row_no)

        # entry = self.edit_entry
        entry = ui.ttk.Entry(self.tree)
        # self.current_entry = entry
        entry.insert(0, content)
        entry.focus_force()
        entry.selection_range(0, 'end')

        def apply_change(event):
            column = self.columns[column_no]
            value = entry.get()
            try:
                column.setvalue(self.atoms, row_no, value)
            except Exception as ex:
                pass
            else:
                text = column.format_value(column.getvalue(self.atoms, row_no))
                self.tree.set(row_id, column_id, value=text)
                self.gui.set_frame()
            finally:
                # pass
                #entry.delete(0, 'end')
                # entry.set('hello')
                self.tree.focus_force()
                entry.destroy()

        entry.bind('<FocusOut>', apply_change)
        entry.bind('<Return>', apply_change)
        entry.bind('<Escape>', lambda *args: entry.destroy())

        x, y, width, height = self.tree.bbox(row_id, column_id)
        # entry.pack()
        entry.place(x=x, y=y, height=height)#width=width)#, height=height)
