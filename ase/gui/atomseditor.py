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
    # We subscribe to gui.draw() calls in order to track changes,
    # but we should have an actual "atoms changed" event instead.

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
                return  # Display error?
            atoms.symbols[i] = value

        self.gui = gui
        self.treeview = treeview
        self._current_entry = None

        columns = []
        symbols_column = Column(
            'symbol', _('symbol'), 60, get_symbol, set_symbol
        )
        columns.append(symbols_column)

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

        for c, axisname in enumerate('xyz'):
            column = Column(
                axisname,
                axisname,
                92,
                GetSetPos(c).get_position,
                GetSetPos(c).set_position,
                format_value=lambda val: f'{val:.4f}',
            )
            columns.append(column)

        self.columns = columns

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
        # Here it is important to leave edit mode since scrolling
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
        self.treeview.delete(*self.treeview.get_children())
        for i in range(len(self.atoms)):
            values = self.get_row_values(i)
            self.treeview.insert(
                '', 'end', text=i, values=values, iid=self.rowid(i)
            )

        selection = np.arange(len(self.atoms))[self.gui.images.selected]

        rowids = [self.rowid(index) for index in selection]
        # Note: selection_set() does *not* fire an event, and therefore
        # we do not need to worry about infinite recursion.
        # However the event listening is wonky now because we need
        # better GUI change listeners.
        self.treeview.selection_set(*rowids)

    def get_row_values(self, i):
        return [
            column.format_value(column.getvalue(self.atoms, i))
            for column in self.columns
        ]

    def define_columns_on_widget(self):
        self.treeview['columns'] = [column.name for column in self.columns]
        for column in self.columns:
            self.treeview.heading(column.name, text=column.displayname)
            self.treeview.column(
                column.name,
                width=column.widget_width,
                anchor='e',
            )

    def rowid(self, rownumber: int) -> str:
        return f'R{rownumber}'

    def rownumber(self, rowid: str) -> int:
        assert rowid.startswith('R'), repr(rowid)
        return int(rowid[1:])

    def set_value(self, column_no: int, row_no: int, value: object) -> None:
        column = self.columns[column_no]
        column.setvalue(self.atoms, row_no, value)
        text = column.format_value(column.getvalue(self.atoms, row_no))

        # The text that we set here is not what matters: It may be rounded.
        # It was column.setvalue() which did the actual change.
        self.treeview.set(self.rowid(row_no), column.name, value=text)

        # (Maybe it is not always necessary to redraw everything.)
        self.gui.set_frame()

    def doubleclick(self, event):
        row_id = self.treeview.identify_row(event.y)
        column_id = self.treeview.identify_column(event.x)
        if not row_id or not column_id:
            return  # clicked outside actual rows/columns
        self.edit_field(row_id, column_id)

    def edit_field(self, row_id, column_id):
        assert column_id.startswith('#'), repr(column_id)
        column_no = int(column_id[1:]) - 1

        if column_no == -1:
            return  # This is the ID column.

        row_no = self.rownumber(row_id)
        assert 0 <= column_no < len(self.columns)
        assert 0 <= row_no < len(self.atoms)

        content = self.columns[column_no].getvalue(self.atoms, row_no)

        assert self._current_entry is None
        entry = ui.ttk.Entry(self.treeview)
        entry.insert(0, content)
        entry.focus_force()
        entry.selection_range(0, 'end')

        def apply_change(_event=None):
            value = entry.get()
            try:
                self.set_value(column_no, row_no, value)
            finally:
                # Focus was given to the text field, now return it:
                self.treeview.focus_force()
                self.leave_edit_mode()

        entry.bind('<FocusOut>', apply_change)
        ui.bind_enter(entry, apply_change)
        entry.bind('<Escape>', lambda *args: self.leave_edit_mode())

        bbox = self.treeview.bbox(row_id, column_id)
        if bbox:  # (bbox is '' when testing without display)
            x, y, width, height = bbox
            entry.place(x=x, y=y, height=height)
        self._current_entry = entry
        return entry, apply_change
