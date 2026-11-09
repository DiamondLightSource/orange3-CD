import numpy as np
import pandas as pd
from AnyQt.QtCore import (
    QAbstractTableModel,
    QModelIndex,
    Qt,
)
from AnyQt.QtWidgets import (
    QFileDialog,
    QHeaderView,
)
from Orange.data import Table
from Orange.widgets import gui
from Orange.widgets.utils.tableview import TableView
from Orange.widgets.widget import Msg, MultiInput, OWWidget

from ...units import Q_
from .utils import DELTA_EPSILON_UNIT, WAVELENGTH_UNIT_KEY


def unit_text(unit):
    """Short display text for a unit string, e.g. ``"µM"``; empty if none."""
    if not unit:
        return ""
    try:
        units = Q_(1, unit).units
        if units == Q_(1, DELTA_EPSILON_UNIT).units:
            return "M\u207b\u00b9 cm\u207b\u00b9"  # the conventional form
        return f"{units:~}"
    except Exception:  # not a unit pint understands: show it as written
        return str(unit)


def column_header(variable, wavelength_unit=None):
    """Column name with the unit of the variable, e.g. ``"Conc [B] (M)"``.

    Spectra tables name their attributes by wavelength, which has no unit of
    its own; it takes the table's wavelength unit instead.
    """
    unit = variable.attributes.get("unit")
    if not unit and wavelength_unit and variable.is_continuous:
        try:
            float(variable.name)
        except ValueError:
            pass
        else:
            unit = wavelength_unit
    text = unit_text(unit)
    return f"{variable.name} ({text})" if text else variable.name


def table_to_dataframe(table):

    columns = []

    blocks = []

    wavelength_unit = table.attributes.get(WAVELENGTH_UNIT_KEY)

    if len(table.domain.attributes):
        columns.extend(
            column_header(var, wavelength_unit) for var in table.domain.attributes
        )

        blocks.append(table.X)

    if len(table.domain.class_vars):
        columns.extend(column_header(var) for var in table.domain.class_vars)

        blocks.append(table.Y.reshape(len(table), -1))

    if len(table.domain.metas):
        columns.extend(column_header(var) for var in table.domain.metas)

        blocks.append(table.metas)

    if not blocks:
        return pd.DataFrame(index=range(len(table)))

    matrix = np.hstack(blocks)

    return pd.DataFrame(
        matrix,
        columns=columns,
    )


METADATA_SHEET = "Metadata"


def metadata_dataframe(sheet_names, tables):
    """One row per table attribute, naming the sheet it belongs to.

    Attributes carry the concentration, pathlength, units and other values
    that the column headers cannot hold. Returns None if there are none.
    """
    rows = [
        (sheet_name, str(key), str(value))
        for sheet_name, table in zip(sheet_names, tables)
        for key, value in table.attributes.items()
    ]
    if not rows:
        return None
    return pd.DataFrame(rows, columns=["Sheet", "Attribute", "Value"])


def make_unique_sheet_names(names):
    """
    Excel sheet names:
      - max 31 chars
      - must be unique
    """

    used = set()
    result = []

    for name in names:
        if not name:
            name = "Sheet"

        # sanitise Excel-invalid characters
        name = (
            str(name)
            .replace("/", "_")
            .replace("\\", "_")
            .replace("*", "_")
            .replace("?", "_")
            .replace("[", "_")
            .replace("]", "_")
            .replace(":", "_")
        ).strip()

        if not name:
            name = "Sheet"

        base = name[:31]
        candidate = base
        counter = 1

        while candidate.lower() in used:
            suffix = f" ({counter})"

            candidate = base[: 31 - len(suffix)] + suffix

            counter += 1

        used.add(candidate.lower())
        result.append(candidate)

    return result


# ----------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------

# Shared invalid index: the default parent of the model's row/column queries.
NO_PARENT = QModelIndex()


class SheetListModel(QAbstractTableModel):
    HEADERS = ("Input", "Sheet Name")

    def __init__(self, parent=None):
        super().__init__(parent)
        self.items = []

    def rowCount(self, parent=NO_PARENT):
        return len(self.items)

    def columnCount(self, parent=NO_PARENT):
        return 2

    def headerData(self, section, orientation, role):

        if role != Qt.DisplayRole:
            return None

        if orientation == Qt.Horizontal:
            return self.HEADERS[section]

        return str(section + 1)

    def data(self, index, role):

        if not index.isValid():
            return None

        row = index.row()
        col = index.column()

        item = self.items[row]

        if role in (Qt.DisplayRole, Qt.EditRole):
            if col == 0:
                return str(row + 1)

            if col == 1:
                return item["sheet_name"]

        return None

    def flags(self, index):

        flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable

        if index.column() == 1:
            flags |= Qt.ItemIsEditable

        return flags

    def setData(self, index, value, role):

        if role == Qt.EditRole and index.column() == 1:
            self.items[index.row()]["sheet_name"] = str(value)

            self.dataChanged.emit(index, index)

            return True

        return False

    def insert_table(self, position, table):

        default_name = getattr(table, "name", None) or f"Sheet{position + 1}"

        self.beginInsertRows(
            QModelIndex(),
            position,
            position,
        )

        self.items.insert(
            position,
            {
                "table": table,
                "sheet_name": default_name,
            },
        )

        self.endInsertRows()

    def update_table(self, position, table):

        if position >= len(self.items):
            self.insert_table(position, table)
            return

        self.items[position]["table"] = table

        self.dataChanged.emit(
            self.index(position, 0),
            self.index(position, 1),
        )

    def remove_table(self, position):

        if position >= len(self.items):
            return

        self.beginRemoveRows(
            QModelIndex(),
            position,
            position,
        )

        del self.items[position]

        self.endRemoveRows()


# ----------------------------------------------------------------------
# Widget
# ----------------------------------------------------------------------


class OWMultiSave(OWWidget):
    name = "Save Workbook"
    description = "Save multiple tables to an Excel workbook in different sheets"
    icon = "icons/SaveWorkbook.svg"

    want_main_area = False

    class Inputs:
        data = MultiInput("Data", Table)

    class Error(OWWidget.Error):
        save_failed = Msg("{}")

    def __init__(self):
        super().__init__()

        self.model = SheetListModel(self)

        box = gui.widgetBox(
            self.controlArea,
            "Workbook Sheets",
        )

        self.view = TableView(self)
        self.view.setModel(self.model)

        header = self.view.horizontalHeader()
        header.setSectionResizeMode(
            0,
            QHeaderView.ResizeToContents,
        )
        header.setSectionResizeMode(
            1,
            QHeaderView.Stretch,
        )

        box.layout().addWidget(self.view)

        gui.button(
            self.controlArea,
            self,
            "Save Workbook",
            callback=self.save_workbook,
        )

    # --------------------------------------------------
    # Inputs
    # --------------------------------------------------

    @Inputs.data.insert
    def insert_data(self, index, data):

        self.model.insert_table(
            index,
            data,
        )

    @Inputs.data
    def set_data(self, index, data):

        self.model.update_table(
            index,
            data,
        )

    @Inputs.data.remove
    def remove_data(self, index):
        self.model.remove_table(index)

    # --------------------------------------------------
    # Save
    # --------------------------------------------------

    def save_workbook(self):
        self.Error.save_failed.clear()

        valid_items = [item for item in self.model.items if item["table"] is not None]
        if not valid_items:
            self.Error.save_failed("There are no tables to save.")
            return

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Workbook",
            "",
            "Excel Workbook (*.xlsx)",
        )

        if not filename:
            return

        if not filename.endswith(".xlsx"):
            filename += ".xlsx"

        tables = [item["table"] for item in valid_items]
        # The metadata sheet is named together with the others so that it can
        # never clash with a table's sheet.
        *sheet_names, metadata_name = make_unique_sheet_names(
            [*(item["sheet_name"] for item in valid_items), METADATA_SHEET]
        )
        metadata = metadata_dataframe(sheet_names, tables)

        try:
            with pd.ExcelWriter(
                filename,
                engine="openpyxl",
            ) as writer:
                for table, sheet_name in zip(tables, sheet_names):
                    table_to_dataframe(table).to_excel(
                        writer,
                        sheet_name=sheet_name,
                        index=False,
                    )
                if metadata is not None:
                    metadata.to_excel(writer, sheet_name=metadata_name, index=False)
        except (OSError, ValueError) as exc:
            self.Error.save_failed(f"Could not save the workbook: {exc}")
