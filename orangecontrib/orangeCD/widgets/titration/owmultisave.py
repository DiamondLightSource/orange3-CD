from pathlib import Path
import numpy as np
import pandas as pd

from AnyQt.QtCore import (
    Qt,
    QModelIndex,
    QAbstractTableModel,
)

from AnyQt.QtWidgets import (
    QFileDialog,
    QHeaderView,
)

from Orange.data import Table
from Orange.widgets.widget import OWWidget, MultiInput
from Orange.widgets import gui
from Orange.widgets.utils.tableview import TableView


def table_to_dataframe(table):

    columns = []

    blocks = []

    if len(table.domain.attributes):

        columns.extend(
            var.name
            for var in table.domain.attributes
        )

        blocks.append(table.X)

    if len(table.domain.class_vars):

        columns.extend(
            var.name
            for var in table.domain.class_vars
        )

        blocks.append(table.Y.reshape(len(table), -1))

    if len(table.domain.metas):

        columns.extend(
            var.name
            for var in table.domain.metas
        )

        blocks.append(table.metas)

    matrix = np.hstack(blocks)

    return pd.DataFrame(
        matrix,
        columns=columns,
    )

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

            candidate = (
                base[: 31 - len(suffix)]
                + suffix
            )

            counter += 1

        used.add(candidate.lower())
        result.append(candidate)

    return result

# ----------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------

class SheetListModel(QAbstractTableModel):

    HEADERS = ["Input", "Sheet Name"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.items = []

    def rowCount(self, parent=QModelIndex()):
        return len(self.items)

    def columnCount(self, parent=QModelIndex()):
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

        if (
            role == Qt.EditRole
            and index.column() == 1
        ):
            self.items[index.row()]["sheet_name"] = str(value)

            self.dataChanged.emit(index, index)

            return True

        return False

    def insert_table(self, position, table):

        default_name = (
            getattr(table, "name", None)
            or f"Sheet{position + 1}"
        )

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

            print(
                "MODEL update -> extending model"
            )

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
    description = (
        "Save multiple tables to an Excel workbook in different sheets"
    )

    want_main_area = False

    class Inputs:
        data = MultiInput("Data", Table)

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

        valid_items = [item for item in self.model.items if item["table"] is not None]
        sheet_names = make_unique_sheet_names([item["sheet_name"] for item in valid_items])

        with pd.ExcelWriter(
            filename,
            engine="openpyxl",
        ) as writer:
            for item, sheet_name in zip(valid_items, sheet_names):
                table = item["table"]

                try:
                    df = table_to_dataframe(table)
                    df.to_excel(
                        writer,
                        sheet_name=sheet_name,
                        index=False,
                    )

                except Exception as ex:
                    raise

