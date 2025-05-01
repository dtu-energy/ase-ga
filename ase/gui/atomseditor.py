from dataclasses import dataclass
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
    # * We need to know whenever anything changes, so we can update the table.
    # * Therefore, the main GUI should provide a listener for any atoms change.
    # * Currently anything that changes the atoms also remembers to call
    #   one or more of draw(), set_frame(), or new_atoms().

    def __init__(self, gui):
        gui.register_vulnerable(self)

        win = ui.Window(_('Edit atoms'), wmtype='utility')

        tree = ui.ttk.Treeview(win.win, selectmode='extended')
        edit_entry = ui.ttk.Entry(win.win)
        edit_entry.pack(side='bottom', fill='x')
        tree.pack(side='left', fill='y')
        bar = ui.ttk.Scrollbar(
            win.win, orient='vertical', command=self.scroll_via_scrollbar
        )
        tree.configure(yscrollcommand=self.scroll_via_treeview)

        tree.column('#0', width=40)
        tree.heading('#0', text=_('id'))

        bar.pack(side='right', fill='y')
        self.scrollbar = bar

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
        self._current_entry = None

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
                    try:
                        value = float(value)
                    except ValueError:
                        return
                    atoms.positions[i, self.c] = value

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

        tree.bind('<Double-1>', self.doubleclick)

        self.edit_entry = edit_entry

    def scroll_via_scrollbar(self, *args, **kwargs):
        self.leave_edit_mode()
        return self.tree.yview(*args, **kwargs)

    def scroll_via_treeview(self, *args, **kwargs):
        # Here it is important to leave edit mode, since scrolling
        # invalidates the widget location.  Alternatively we could keep
        # it open as long as we move it but that sounds like work
        self.leave_edit_mode()
        return self.scrollbar.set(*args, **kwargs)

    def leave_edit_mode(self):
        if self._current_entry is not None:
            self._current_entry.destroy()
            self._current_entry = None
            self.tree.focus_force()

    @property
    def atoms(self):
        return self.gui.atoms

    def update_table_from_atoms(self):
        for i in range(len(self.atoms)):
            values = self.get_row_values(i)
            self.tree.insert('', 'end', text=i, values=values)
        self.add_columns_to_widget()

    def get_row_values(self, i):
        return [
            column.format_value(column.getvalue(self.gui.atoms, i))
            for column in self.columns
        ]

    def define_column(self, *args, **kwargs):
        column = Column(*args, **kwargs)
        self.columns.append(column)

    def add_columns_to_widget(self):
        self.tree['columns'] = [column.name for column in self.columns]
        for column in self.columns:
            self.tree.heading(column.name, text=column.displayname)
            self.tree.column(column.name, width=column.widget_width, anchor='e')

    def doubleclick(self, event):
        row_id = self.tree.identify_row(event.y)
        column_id = self.tree.identify_column(event.x)

        assert column_id.startswith('#'), repr(column_id)
        assert row_id.startswith('I'), repr(row_id)

        column_no = int(column_id[1:]) - 1
        if column_no == -1:
            return  # This is the ID column.

        # WTH why in the name of the devil does it use hexadecimal
        row_no = int(row_id[1:], base=16) - 1

        assert 0 <= column_no < len(self.columns)
        assert 0 <= row_no < len(self.atoms)

        content = self.columns[column_no].getvalue(self.atoms, row_no)

        assert self._current_entry is None
        entry = ui.ttk.Entry(self.tree)
        entry.insert(0, content)
        entry.focus_force()
        entry.selection_range(0, 'end')

        def apply_change(event):
            column = self.columns[column_no]
            value = entry.get()
            try:
                column.setvalue(self.atoms, row_no, value)
                text = column.format_value(column.getvalue(self.atoms, row_no))
                self.tree.set(row_id, column_id, value=text)
                self.gui.set_frame()
            finally:
                self.tree.focus_force()
                self.leave_edit_mode()

        entry.bind('<FocusOut>', apply_change)
        entry.bind('<Return>', apply_change)
        entry.bind('<Escape>', lambda *args: self.leave_edit_mode())

        x, y, width, height = self.tree.bbox(row_id, column_id)
        entry.place(x=x, y=y, height=height)
        self._current_entry = entry
