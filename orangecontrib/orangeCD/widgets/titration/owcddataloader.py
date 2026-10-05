"""Load circular-dichroism CSV files into a spectra table for Orange."""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
import pyqtgraph as pg
from AnyQt.QtCore import Qt
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Msg, Output, OWWidget

from .utils import (
    SPECTRUM_META,
    SPECTRUM_UNIT_KEY,
    WAVELENGTH_UNIT_KEY,
    unit_string,
    unit_symbol,
)

WAVELENGTH_UNITS = ["nanometer", "angstrom", "micrometer"]
CD_SIGNAL_UNITS = ["millidegree", "degree"]


def test_empty_line(line: str) -> str | None:
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
                key = tokens[0].split("#", maxsplit=1)[1].strip()
            except IndexError:
                comments[line.strip("#\n")] = None
            else:
                comments[key] = tokens[1].strip().split(",")[0]
        else:
            comments[line.strip("#\n")] = None
    return comments


def parse_data(lines: list[str]) -> dict[str, pd.DataFrame]:
    parsed: dict[str, pd.DataFrame] = {}
    sections = [i for i, line in enumerate(lines) if "Wavelength" in line]
    for start, stop in zip(sections, sections[1:]):
        title = lines[start - 1].split(",")[0].strip()
        rows = [
            [float(value) for value in line.strip().split(",") if value]
            for line in lines[start + 2 : stop - 2]
            if line.strip(",\n ")
        ]
        if not rows:
            continue
        values = np.asarray(rows, dtype=float)
        if values.ndim != 2 or values.shape[1] < 2:
            raise ValueError(f"Invalid data in section {title!r}")
        measurements = values[:, 1:]
        output = np.column_stack((measurements, measurements.mean(axis=1)))
        frame = pd.DataFrame(
            output,
            index=values[:, 0],
            columns=[str(i) for i in range(measurements.shape[1])] + ["Average"],
        )
        frame.index.name = "Wavelength"
        parsed[title] = frame
    return parsed


def file_parser(filename: str) -> dict[str, object]:
    with open(filename, encoding="utf-8") as stream:
        lines = stream.readlines()
    starts = [
        i for i, line in enumerate(lines)
        if ":" in line.split(",")[0] and "#" not in line.split(",")[0]
    ]
    starts.append(len(lines))
    parsers = {"Remarks": parse_remarks, "Data": parse_data}
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
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", Path(path).name)
    ]


def dataframe_to_orange_table(
    dataframe: pd.DataFrame, wavelength_unit: str, spectrum_unit: str
) -> Table:
    """Convert a wavelength-indexed frame into a spectra-as-rows table.

    Columns of ``dataframe`` (one per spectrum) become rows; the index
    (wavelength) becomes the attributes, in ascending order.
    """
    output = dataframe.copy()
    output.columns = [
        " | ".join(map(str, col)) if isinstance(col, tuple) else str(col)
        for col in output.columns
    ]
    # Ascending wavelength axis, whatever order the instrument wrote.
    output = output.sort_index()

    domain = Domain(
        [ContinuousVariable(str(wavelength)) for wavelength in output.index.values],
        metas=[StringVariable(SPECTRUM_META)],
    )
    table = Table.from_numpy(
        domain,
        output.to_numpy(dtype=float).T,
        metas=np.array(output.columns, dtype=object)[:, np.newaxis],
    )
    table.attributes[WAVELENGTH_UNIT_KEY] = unit_string(wavelength_unit)
    table.attributes[SPECTRUM_UNIT_KEY] = unit_string(spectrum_unit)
    return table


class OWCDDataLoader(OWWidget):
    name = "CD Data Loader"
    description = (
        "Load CD spectra from CSV files into a table of spectra, ready for "
        "spectral preprocessing."
    )
    icon = "icons/CDDataLoader.svg"
    priority = 20
    want_main_area = True
    resizing_enabled = True

    solution_a_file = Setting("")
    solution_b_file = Setting("")
    buffer_file = Setting("")
    data_files = Setting([])
    selected_data_files = Setting([])
    wavelength_unit = Setting(WAVELENGTH_UNITS[0])
    cd_unit = Setting(CD_SIGNAL_UNITS[0])
    data_file_names: list[str] = []

    class Outputs:
        data = Output("CD Data", Table)

    class Error(OWWidget.Error):
        load_failed = Msg("{}")

    def __init__(self) -> None:
        super().__init__()
        self.cd_data: pd.DataFrame | None = None
        self._build_controls()
        self._build_plot()
        self._refresh_file_list()
        if self.data_files:
            self._load_and_plot_raw_data()
        self.commit()

    def _build_controls(self) -> None:
        references = gui.widgetBox(self.controlArea, "Reference files (optional)")
        for label, setting in (
            ("Solution A", "solution_a_file"),
            ("Solution B", "solution_b_file"),
            ("Buffer", "buffer_file"),
        ):
            self._add_reference_selector(references, label, setting)

        data_box = gui.widgetBox(self.controlArea, "Titration data files")
        self.data_file_list = gui.listBox(
            data_box,
            self,
            "selected_data_files",
            "data_file_names",
            selectionMode=gui.QtWidgets.QAbstractItemView.ExtendedSelection,
        )
        self.data_file_list.setMinimumWidth(480)
        self.data_file_list.setMinimumHeight(150)
        buttons = gui.hBox(data_box)
        gui.button(buttons, self, "Select files…", callback=self._choose_data_files)
        gui.button(buttons, self, "Clear", callback=self._clear_data_files)

        units_box = gui.widgetBox(self.controlArea, "Units of the data files")
        gui.comboBox(
            units_box, self, "wavelength_unit", label="Wavelength:",
            items=WAVELENGTH_UNITS, sendSelectedValue=True,
            orientation=Qt.Horizontal, callback=self._units_changed,
        )
        gui.comboBox(
            units_box, self, "cd_unit", label="CD signal:",
            items=CD_SIGNAL_UNITS, sendSelectedValue=True,
            orientation=Qt.Horizontal, callback=self._units_changed,
        )
        gui.rubber(self.controlArea)

    def _apply_unit_labels(self) -> None:
        self.plot_widget.setLabel(
            "bottom", "Wavelength", units=unit_symbol(self.wavelength_unit)
        )
        self.plot_widget.setLabel(
            "left", "Circular dichroism", units=unit_symbol(self.cd_unit)
        )

    def _units_changed(self) -> None:
        self._apply_unit_labels()
        self.commit()

    def _add_reference_selector(self, parent, label: str, setting: str) -> None:
        row = gui.hBox(parent)
        editor = gui.lineEdit(
            row, self, setting, label=label, orientation=Qt.Horizontal,
            callback=self.commit,
        )
        editor.setMinimumWidth(380)
        editor.setToolTip(str(getattr(self, setting)))
        editor.textChanged.connect(editor.setToolTip)
        gui.button(
            row, self, "Browse…",
            callback=lambda _checked=False, name=setting: self._choose_reference(name),
        )

    def _build_plot(self) -> None:
        plot_box = gui.vBox(self.mainArea)
        gui.widgetLabel(plot_box, "Loaded spectra")
        self.plot_widget = pg.PlotWidget(plot_box)
        self.plot_widget.showGrid(x=True, y=True, alpha=0.2)
        self._apply_unit_labels()
        plot_box.layout().addWidget(self.plot_widget)

    def _choose_reference(self, setting: str) -> None:
        filename, _ = gui.QtWidgets.QFileDialog.getOpenFileName(
            self, "Select reference CSV", "", "CSV files (*.csv);;All files (*)"
        )
        if filename:
            setattr(self, setting, filename)
            self.commit()

    def _choose_data_files(self) -> None:
        filenames, _ = gui.QtWidgets.QFileDialog.getOpenFileNames(
            self, "Select titration CSV files", "", "CSV files (*.csv);;All files (*)"
        )
        if filenames:
            self.data_files = sorted(dict.fromkeys(filenames), key=natural_key)
            self._refresh_file_list()
            self._load_and_plot_raw_data()
            self.commit()

    def _clear_data_files(self) -> None:
        self.data_files = []
        self.cd_data = None
        self._refresh_file_list()
        self._refresh_plot()
        self.commit()

    def _refresh_file_list(self) -> None:
        self.data_file_names = [Path(path).name for path in self.data_files]
        self.selected_data_files = list(range(len(self.data_file_names)))

    def _load_and_plot_raw_data(self) -> None:
        self.Error.load_failed.clear()
        try:
            self.cd_data = pd.concat(
                [average_cd_series(path) for path in self.data_files], axis=1
            )
        except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
            self.cd_data = None
            self.Error.load_failed(str(exc))
        self._refresh_plot()

    def _refresh_plot(self) -> None:
        self.plot_widget.clear()
        if self.cd_data is not None:
            x = self.cd_data.index.to_numpy(dtype=float)
            count = max(len(self.cd_data.columns), 1)
            for index, name in enumerate(self.cd_data.columns):
                self.plot_widget.plot(
                    x, self.cd_data[name].to_numpy(dtype=float),
                    pen=pg.mkPen(pg.intColor(index, count), width=1.2),
                )
        self.plot_widget.enableAutoRange()

    def commit(self) -> None:
        """Send every loaded spectrum, uncorrected, as one table.

        Reference files become ``Background | <role>`` rows and data files
        ``<file name> | raw_data`` rows.
        """
        self.Error.load_failed.clear()
        references = {
            "sol_A": self.solution_a_file,
            "sol_B": self.solution_b_file,
            "buffer": self.buffer_file,
        }
        try:
            series: dict[tuple[str, str], pd.Series] = {}
            for role, path in references.items():
                if path:
                    series[("Background", role)] = average_cd_series(path)
            for path in self.data_files:
                key = (Path(path).name, "raw_data")
                if key in series:
                    raise ValueError(f"Duplicate file name: {key[0]}")
                series[key] = average_cd_series(path)
            if not series:
                self.Outputs.data.send(None)
                return
            frame = pd.DataFrame(series)
            if frame.isna().any().any():
                raise ValueError(
                    "The selected files do not all share the same wavelength axis"
                )
            table = dataframe_to_orange_table(
                frame, self.wavelength_unit, self.cd_unit
            )
        except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
            self.Error.load_failed(str(exc))
            self.Outputs.data.send(None)
            return
        self.Outputs.data.send(table)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWCDDataLoader).run()
