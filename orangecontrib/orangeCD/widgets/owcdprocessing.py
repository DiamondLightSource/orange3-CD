"""Orange widget for reading and processing circular-dichroism titration CSVs."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from AnyQt.QtCore import QAbstractTableModel, QModelIndex, Qt
from AnyQt.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QTableView,
    QVBoxLayout,
    QWidget,
)
from Orange.data import ContinuousVariable, Domain, Table
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget


def test_empty_line(line: str) -> str | None:
    """Return cleaned non-empty comma-separated text, or None."""
    if not "".join(line.strip().split(",")):
        return None
    return ", ".join(line.strip().split(","))


def parse_remarks(lines: Iterable[str]) -> dict[str, str | None]:
    comments: dict[str, str | None] = {}
    for line in lines:
        tokens = line.split(":")
        if len(tokens) == 1:
            cleaned = test_empty_line(line)
            if cleaned is not None:
                comments["Header"] = cleaned
        elif len(tokens) == 2:
            try:
                entry_key = tokens[0].split("#", maxsplit=1)[1].strip()
            except IndexError:
                comments[line.strip("#\n")] = None
            else:
                comments[entry_key] = tokens[1].strip().split(",")[0]
        else:
            comments[line.strip("#\n")] = None
    return comments


def parse_data(lines: list[str]) -> dict[str, pd.DataFrame]:
    """Parse all spectral data sections from a CD CSV file."""
    data_dict: dict[str, pd.DataFrame] = {}
    sections = [index for index, line in enumerate(lines) if "Wavelength" in line]

    for start, stop in zip(sections, sections[1:]):
        section_title = lines[start - 1].split(",")[0].strip()
        section_lines = lines[start + 2 : stop - 2]
        rows = [
            [float(value) for value in line.strip().split(",") if value]
            for line in section_lines
            if line.strip(",\n ")
        ]
        if not rows:
            continue

        values = np.asarray(rows, dtype=float)
        if values.ndim != 2 or values.shape[1] < 2:
            raise ValueError(f"Invalid data in section {section_title!r}")

        measurements = values[:, 1:]
        output = np.column_stack((measurements, measurements.mean(axis=1)))
        dataframe = pd.DataFrame(
            output,
            index=values[:, 0],
            columns=[str(i) for i in range(measurements.shape[1])] + ["Average"],
        )
        dataframe.index.name = "Wavelength"
        data_dict[section_title] = dataframe

    return data_dict


def file_parser(filename: str) -> dict[str, object]:
    """Parse one instrument CSV file into its named sections."""
    with open(filename, encoding="utf-8") as stream:
        lines = stream.readlines()

    starts = [
        index
        for index, line in enumerate(lines)
        if ":" in line.split(",")[0] and "#" not in line.split(",")[0]
    ]
    starts.append(len(lines))
    """
    TODO: add more parsers for sections of the input data 
    """
    parsers = {"Remarks": parse_remarks, 
               "Data": parse_data}
    parsed: dict[str, object] = {}

    for start, stop in zip(starts, starts[1:]):
        title = lines[start].split(":", maxsplit=1)[0]
        parser = parsers.get(title)
        if parser is not None:
            parsed[title] = parser(lines[start + 1 : stop])

    return parsed


def average_cd_series(filename: str) -> pd.Series:
    parsed = file_parser(filename)
    try:
        series = parsed["Data"]["CircularDichroism"]["Average"]
    except (KeyError, TypeError) as exc:
        raise ValueError(
            f"{Path(filename).name}: CircularDichroism/Average was not found"
        ) from exc
    return series.rename(Path(filename).name)


def natural_key(path: str) -> list[object]:
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", Path(path).name)]


def orange_table_to_dataframe(table: Table) -> pd.DataFrame:
    """Convert an Orange table's variables and metas to a pandas DataFrame."""
    columns: dict[str, np.ndarray] = {}
    for variable in (*table.domain.attributes, *table.domain.class_vars,
                     *table.domain.metas):
        columns[variable.name] = table.get_column(variable)
    return pd.DataFrame(columns)


def flattened_dataframe(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Return a display-ready copy with flattened column names."""
    output = dataframe.copy()
    output.columns = [
        " | ".join(map(str, column)) if isinstance(column, tuple) else str(column)
        for column in output.columns
    ]
    output.index.name = "Wavelength"
    return output.reset_index()


def dataframe_to_orange_table(dataframe: pd.DataFrame) -> Table:
    """Create an Orange table with wavelength stored as a continuous meta."""
    output = flattened_dataframe(dataframe)
    wavelength = ContinuousVariable("Wavelength")
    attribute_names = [column for column in output.columns if column != "Wavelength"]
    domain = Domain(
        [ContinuousVariable(column) for column in attribute_names],
        metas=[wavelength],
    )
    x = output[attribute_names].to_numpy(dtype=float)
    metas = output[["Wavelength"]].to_numpy(dtype=float)
    return Table.from_numpy(domain, x, metas=metas)


class DataFrameModel(QAbstractTableModel):
    """Read-only model used for the processed-data preview."""

    def __init__(self) -> None:
        super().__init__()
        self._dataframe = pd.DataFrame()

    def set_dataframe(self, dataframe: pd.DataFrame) -> None:
        self.beginResetModel()
        self._dataframe = dataframe.copy()
        self.endResetModel()

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._dataframe.index)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._dataframe.columns)

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole):
        if not index.isValid():
            return None
        value = self._dataframe.iat[index.row(), index.column()]
        if role == Qt.DisplayRole:
            if pd.isna(value):
                return ""
            if isinstance(value, (float, np.floating)):
                return f"{value:g}"
            return str(value)
        if role == Qt.TextAlignmentRole and isinstance(
            value, (int, float, np.integer, np.floating)
        ):
            return int(Qt.AlignRight | Qt.AlignVCenter)
        return None

    def headerData(
        self, section: int, orientation: Qt.Orientation, role: int = Qt.DisplayRole
    ):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            return str(self._dataframe.columns[section])
        return str(section + 1)


class OWCDTitrationProcessing(OWWidget):
    name = "CD Titration Processing"
    description = "Read CD spectra and apply titration reference corrections."
    icon = "icons/Titration.svg"
    priority = 20
    want_main_area = False
    resizing_enabled = True

    solution_a_file = Setting("")
    solution_b_file = Setting("")
    buffer_file = Setting("")
    data_files = Setting([])

    class Inputs:
        titration = Input("Titration Table", Table)

    class Outputs:
        data = Output("Processed CD Data", Table)

    class Error(OWWidget.Error):
        invalid_input = Msg("{}")
        processing_failed = Msg("{}")

    class Warning(OWWidget.Warning):
        no_titration = Msg("Connect a Titration Table before processing.")

    def __init__(self) -> None:
        super().__init__()
        self.titration: Table | None = None
        self._path_edits: dict[str, QLineEdit] = {}
        self._build_controls()
        # self._build_main_area()
        self._refresh_file_list()

    def _build_controls(self) -> None:
        references = QGroupBox("Reference files", self.controlArea)
        reference_form = QFormLayout(references)
        self.controlArea.layout().addWidget(references)

        for label, setting_name in (
            ("Solution A", "solution_a_file"),
            ("Solution B", "solution_b_file"),
            ("Buffer", "buffer_file"),
        ):
            row = QWidget(references)
            layout = QHBoxLayout(row)
            layout.setContentsMargins(0, 0, 0, 0)
            edit = QLineEdit(str(getattr(self, setting_name)), row)
            edit.setMinimumWidth(380)
            edit.setToolTip(edit.text())
            edit.textChanged.connect(edit.setToolTip)
            edit.editingFinished.connect(self._read_reference_fields)
            button = QPushButton("Browse…", row)
            button.clicked.connect(
                lambda _checked=False, name=setting_name: self._choose_reference(name)
            )
            layout.addWidget(edit, 1)
            layout.addWidget(button)
            reference_form.addRow(label, row)
            self._path_edits[setting_name] = edit

        data_box = QGroupBox("Titration data files", self.controlArea)
        data_layout = QVBoxLayout(data_box)
        self.data_file_list = QListWidget(data_box)
        self.data_file_list.setMinimumHeight(150)
        data_layout.addWidget(self.data_file_list)

        buttons = QHBoxLayout()
        choose = QPushButton("Select files…", data_box)
        choose.clicked.connect(self._choose_data_files)
        clear = QPushButton("Clear", data_box)
        clear.clicked.connect(self._clear_data_files)
        buttons.addWidget(choose)
        buttons.addWidget(clear)
        data_layout.addLayout(buttons)
        self.controlArea.layout().addWidget(data_box)

        process_button = QPushButton("Process data", self.controlArea)
        process_button.clicked.connect(self.process)
        self.controlArea.layout().addWidget(process_button)
        self.controlArea.layout().addStretch(1)

    """
    Initially had this as a kind of viewport to live preview the output data table, but it took up too much processing power

    could be added back in if wanted
    """
    # def _build_main_area(self) -> None:
    #     container = QWidget(self.mainArea)
    #     layout = QVBoxLayout(container)
    #     self.status_label = QLabel("No processed output", container)
    #     self.status_label.setWordWrap(True)

    #     self.preview_model = DataFrameModel()
    #     self.preview = QTableView(container)
    #     self.preview.setModel(self.preview_model)
    #     self.preview.setAlternatingRowColors(True)
    #     self.preview.setEditTriggers(QAbstractItemView.NoEditTriggers)
    #     self.preview.setSelectionBehavior(QAbstractItemView.SelectRows)
    #     self.preview.horizontalHeader().setSectionResizeMode(
    #         QHeaderView.ResizeToContents
    #     )
    #     self.preview.horizontalHeader().setStretchLastSection(True)

    #     layout.addWidget(self.status_label)
    #     layout.addWidget(self.preview, 1)
    #     self.mainArea.layout().addWidget(container)

    def _choose_reference(self, setting_name: str) -> None:
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select reference CSV", "", "CSV files (*.csv);;All files (*)"
        )
        if filename:
            setattr(self, setting_name, filename)
            self._path_edits[setting_name].setText(filename)

    def _choose_data_files(self) -> None:
        filenames, _ = QFileDialog.getOpenFileNames(
            self, "Select titration CSV files", "", "CSV files (*.csv);;All files (*)"
        )
        if filenames:
            self.data_files = sorted(dict.fromkeys(filenames), key=natural_key)
            self._refresh_file_list()

    def _clear_data_files(self) -> None:
        self.data_files = []
        self._refresh_file_list()
        self.Outputs.data.send(None)

    def _refresh_file_list(self) -> None:
        self.data_file_list.clear()
        self.data_file_list.addItems([Path(path).name for path in self.data_files])

    def _read_reference_fields(self) -> None:
        for name, edit in self._path_edits.items():
            setattr(self, name, edit.text().strip())

    @Inputs.titration
    def set_titration(self, table: Table | None) -> None:
        self.titration = table
        if table is not None and self.data_files:
            self.process()
        elif table is None:
            self.Outputs.data.send(None)

    def process(self) -> None:
        self.Error.clear()
        self.Warning.clear()
        self._read_reference_fields()

        if self.titration is None:
            self.Warning.no_titration()
            self.Outputs.data.send(None)
            return

        references = {
            "sol_A": self.solution_a_file,
            "sol_B": self.solution_b_file,
            "buffer": self.buffer_file,
        }
        missing = [name for name, path in references.items() if not path]
        if missing:
            self.Error.invalid_input(
                "Select all reference files; missing: " + ", ".join(missing)
            )
            self.Outputs.data.send(None)
            return
        if not self.data_files:
            self.Error.invalid_input("Select at least one titration data file")
            self.Outputs.data.send(None)
            return

        try:
            titration = orange_table_to_dataframe(self.titration)
            required = {"dilution_factor", "normalised_molar_ratio"}
            absent = required.difference(titration.columns)
            if absent:
                raise ValueError(
                    "Titration input is missing columns: " + ", ".join(sorted(absent))
                )
            if len(titration) != len(self.data_files):
                raise ValueError(
                    f"The titration table has {len(titration)} rows but "
                    f"{len(self.data_files)} data files were selected"
                )

            reference_data = pd.concat(
                {name: average_cd_series(path) for name, path in references.items()},
                axis=1,
            )
            reference_data["sol_A_buffer_subtracted"] = (
                reference_data["sol_A"] - reference_data["buffer"]
            )
            reference_data["sol_B_buffer_subtracted"] = (
                reference_data["sol_B"] - reference_data["buffer"]
            )
            reference_data["sol_A_buffer_subtracted_zeroed"] = (
                reference_data["sol_A_buffer_subtracted"]
                - reference_data["sol_A_buffer_subtracted"].loc[:450].mean()
            )

            cd_data = pd.concat(
                [average_cd_series(path) for path in self.data_files], axis=1
            )
            if cd_data.isna().any().any() or reference_data.isna().any().any():
                raise ValueError(
                    "The selected files do not all share the same wavelength axis"
                )

            processed: dict[tuple[str, str], pd.Series | float] = {
                ("Background", column): reference_data[column].round(3)
                for column in reference_data.columns
            }
            for index, series_name in enumerate(cd_data.columns):
                raw = cd_data[series_name]

                buffer_subtraction = (
                    raw - reference_data["buffer"]
                ) * float(titration.iloc[index]["dilution_factor"])

                sol_a_subtraction = (
                    buffer_subtraction
                    - reference_data["sol_A_buffer_subtracted"]
                )

                subtract_fraction_sol_b = (
                    sol_a_subtraction
                    - (reference_data["sol_B_buffer_subtracted"]
                    * float(titration.iloc[index]["normalised_molar_ratio"])
                        )
                )

                """
                # TODO: 
                #   this should be done properly, not hard coded into some flat region like this.
                #   However, this approximately mimics what CDApps does at the moment:
                #       1) define some region in the flat part of the signal
                #       2) take an average 
                #       3) subtract from across the whole spectrum 
                #   may take some future unpicking, but hopefully not too much...
                """
                zero_level = subtract_fraction_sol_b.loc[:450].mean()

                plus_sol_a = (
                    zero_level
                    + reference_data["sol_A_buffer_subtracted_zeroed"]
                )

                processed.update({
                    (series_name, "raw_data"): raw.round(3),
                    (series_name, "buffer_subtraction"): buffer_subtraction.round(3),
                    (series_name, "sol_A_subtraction"): sol_a_subtraction.round(3),
                    (series_name, "subtract_frac_sol_B"): subtract_fraction_sol_b.round(3),
                    (series_name, "set_to_zero"): round(float(zero_level), 3),
                    (series_name, "plus_sol_A"): plus_sol_a.round(3),
                })

            output = pd.DataFrame(processed, index=cd_data.index)
            orange_output = dataframe_to_orange_table(output)
        except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
            self.Error.processing_failed(str(exc))
            # self.status_label.setText("Processing failed")
            # self.preview_model.set_dataframe(pd.DataFrame())
            self.Outputs.data.send(None)
            return

        # self.status_label.setText(
        #     f"Processed {len(self.data_files)} spectra at "
        #     f"{len(output.index)} wavelengths; output has {len(output.columns) + 1} columns."
        # )
        # self.preview_model.set_dataframe(flattened_dataframe(output))
        self.Outputs.data.send(orange_output)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWCDTitrationProcessing).run()