from dataclasses import dataclass
from typing import Callable

import numpy as np

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
        gui.obs.change_atoms.register(self.update_table_from_atoms)

        win = ui.Window(_('Edit atoms'), wmtype='utility')

        treeview = ui.ttk.Treeview(win.win, selectmode='extended')
        edit_entry = ui.ttk.Entry(win.win)
        edit_entry.pack(side='bottom', fill='x')
        treeview.pack(side='left', fill='y')
        bar = ui.ttk.Scrollbar(
            win.win, orient='vertical', command=self.scroll_via_scrollbar
        )
        treeview.configure(yscrollcommand=self.scroll_via_treeview)

        treeview.column('#0', width=40)
        treeview.heading('#0', text=_('id'))

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
        self.treeview = treeview
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

        treeview.bind('<Double-1>', self.doubleclick)
        treeview.bind('<<TreeviewSelect>>', self.treeview_selection_changed)

        self.define_columns_on_widget()
        self.update_table_from_atoms()

        self.edit_entry = edit_entry

    def treeview_selection_changed(self, event):
        selected_items = self.treeview.selection()
        indices = [self.rownumber(item) for item in selected_items]
        self.gui.set_selected_atoms(indices)

    def scroll_via_scrollbar(self, *args, **kwargs):
        self.leave_edit_mode()
        return self.treeview.yview(*args, **kwargs)

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
            self.treeview.focus_force()

    @property
    def atoms(self):
        return self.gui.atoms

    def update_table_from_atoms(self):
        previous_selection = [
            self.rownumber(rowid) for rowid in self.treeview.selection()
        ]

        self.treeview.delete(*self.treeview.get_children())
        for i in range(len(self.atoms)):
            values = self.get_row_values(i)
            self.treeview.insert(
                '', 'end', text=i, values=values, iid=self.rowid(i)
            )

        selection = np.arange(len(self.atoms))[self.gui.images.selected]

        rowids = [self.rowid(index) for index in selection]
        # Note: This will cause us to update the selection, which can lead
        # to infinite recursion since we'll also end up updating the table.
        # Except the GUI catches that the selection did not in fact change,
        # and therefore does not update again.
        #
        # This could be made a little bit more obvious/robust
        self.treeview.selection_set(*rowids)

    def get_row_values(self, i):
        return [
            column.format_value(column.getvalue(self.gui.atoms, i))
            for column in self.columns
        ]

    def define_column(self, *args, **kwargs):
        column = Column(*args, **kwargs)
        self.columns.append(column)

    def define_columns_on_widget(self):
        self.treeview['columns'] = [column.name for column in self.columns]
        for column in self.columns:
            self.treeview.heading(column.name, text=column.displayname)
            self.treeview.column(
                column.name, width=column.widget_width, anchor='e'
            )

    def rowid(self, rownumber: int) -> str:
        return f'R{rownumber}'

    def rownumber(self, rowid: str) -> int:
        assert rowid[0] == 'R'
        return int(rowid[1:])

    def doubleclick(self, event):
        row_id = self.treeview.identify_row(event.y)
        column_id = self.treeview.identify_column(event.x)

        assert column_id.startswith('#'), repr(column_id)

        column_no = int(column_id[1:]) - 1
        if column_no == -1:
            return  # This is the ID column.

        row_no = self.rownumber(rowid)

        assert 0 <= column_no < len(self.columns)
        assert 0 <= row_no < len(self.atoms)

        content = self.columns[column_no].getvalue(self.atoms, row_no)

        assert self._current_entry is None
        entry = ui.ttk.Entry(self.treeview)
        entry.insert(0, content)
        entry.focus_force()
        entry.selection_range(0, 'end')

        def apply_change(event):
            column = self.columns[column_no]
            value = entry.get()
            try:
                column.setvalue(self.atoms, row_no, value)
                text = column.format_value(column.getvalue(self.atoms, row_no))
                self.treeview.set(row_id, column_id, value=text)
                self.gui.set_frame()
            finally:
                self.treeview.focus_force()
                self.leave_edit_mode()

        entry.bind('<FocusOut>', apply_change)
        entry.bind('<Return>', apply_change)
        entry.bind('<Escape>', lambda *args: self.leave_edit_mode())

        x, y, width, height = self.treeview.bbox(row_id, column_id)
        entry.place(x=x, y=y, height=height)
        self._current_entry = entry
