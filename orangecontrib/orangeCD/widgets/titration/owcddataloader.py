"""Load circular-dichroism CSV files into a spectra table for Orange."""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd
import pyqtgraph as pg
from AnyQt.QtCore import Qt
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from .utils import (
    POINT_META,
    ROLE_BUFFER,
    ROLE_META,
    ROLE_SOL_A,
    ROLE_SOL_B,
    ROLE_TITRATION,
    SAMPLE_META,
    SPECTRUM_META,
    STAGE_META,
    STAGE_RAW,
    SpectrumLabel,
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
    for start, stop in pairwise(sections):
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
        i
        for i, line in enumerate(lines)
        if ":" in line.split(",")[0] and "#" not in line.split(",")[0]
    ]
    starts.append(len(lines))
    parsers = {"Remarks": parse_remarks, "Data": parse_data}
    parsed: dict[str, object] = {}
    for start, stop in pairwise(starts):
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


_NUMBER = re.compile(r"(\d+(?:\.\d+)?)")


def natural_key(path: str) -> list[object]:
    """Sort key that orders numbers in a file name by value.

    Decimals are read as one number, so ``0.25`` sorts before ``0.3``.
    """
    return [
        float(part) if index % 2 else part.lower()
        for index, part in enumerate(_NUMBER.split(Path(path).name))
    ]


def dataframe_to_orange_table(
    dataframe: pd.DataFrame,
    wavelength_unit: str,
    spectrum_unit: str,
    points: Sequence[float] | None = None,
) -> Table:
    """Convert a wavelength-indexed frame into a spectra-as-rows table.

    Columns of ``dataframe`` (one per spectrum) become rows; the index
    (wavelength) becomes the attributes, in ascending order. A column is
    either a ``(sample, stage, role)`` tuple, which gives the ``Spectrum``,
    ``Sample``, ``Stage`` and ``Role`` metas, or a plain name, which gives
    only ``Spectrum``. ``points``, if given, holds the titration point number
    of each column (NaN for none) and becomes the ``Point`` meta.
    """
    output = dataframe.copy()
    labelled = len(output.columns) > 0 and all(
        isinstance(col, tuple) and len(col) == 3 for col in output.columns
    )
    labels = [SpectrumLabel(*col) for col in output.columns] if labelled else []
    output.columns = (
        [label.name for label in labels]
        if labelled
        else [
            " | ".join(map(str, col)) if isinstance(col, tuple) else str(col)
            for col in output.columns
        ]
    )
    # Ascending wavelength axis, whatever order the instrument wrote.
    output = output.sort_index()

    meta_variables = [StringVariable(SPECTRUM_META)]
    meta_columns = [np.array(output.columns, dtype=object)]
    if labelled:
        for meta, attribute in (
            (SAMPLE_META, "sample"),
            (STAGE_META, "stage"),
            (ROLE_META, "role"),
        ):
            meta_variables.append(StringVariable(meta))
            meta_columns.append(
                np.array([getattr(label, attribute) for label in labels], dtype=object)
            )
    if points is not None:
        meta_variables.append(ContinuousVariable(POINT_META))
        meta_columns.append(np.asarray(points, dtype=float).astype(object))
    domain = Domain(
        [ContinuousVariable(str(wavelength)) for wavelength in output.index.values],
        metas=meta_variables,
    )
    table = Table.from_numpy(
        domain,
        output.to_numpy(dtype=float).T,
        metas=np.column_stack(meta_columns),
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

    class Inputs:
        titration = Input("Titration Table", Table)

    class Outputs:
        data = Output("CD Data", Table)

    class Error(OWWidget.Error):
        load_failed = Msg("{}")

    class Warning(OWWidget.Warning):
        point_count_mismatch = Msg(
            "There are {files} titration data files but the Titration Table has "
            "{rows} rows, so the files are not numbered."
        )
        cropped_wavelengths = Msg(
            "The files cover different wavelengths. Only the {count} wavelengths "
            "common to all files ({low:g} to {high:g}) are used."
        )

    def __init__(self) -> None:
        super().__init__()
        # List model displayed by gui.listBox; must exist before the controls.
        self.data_file_names: list[str] = []
        self.cd_data: pd.DataFrame | None = None
        self.titration: Table | None = None
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
        self.data_file_list.setToolTip(
            "Every file in the list is loaded, in this order. Selecting files "
            "only lets you move them with the buttons below."
        )
        self.data_file_list.setMinimumWidth(480)
        self.data_file_list.setMinimumHeight(150)
        buttons = gui.hBox(data_box)
        gui.button(buttons, self, "Select files…", callback=self._choose_data_files)
        gui.button(
            buttons,
            self,
            "Move up",
            callback=lambda: self._move_selected(-1),
            tooltip="Move the selected files up the list.",
        )
        gui.button(
            buttons,
            self,
            "Move down",
            callback=lambda: self._move_selected(1),
            tooltip="Move the selected files down the list.",
        )
        gui.button(buttons, self, "Clear", callback=self._clear_data_files)

        units_box = gui.widgetBox(self.controlArea, "Units of the data files")
        gui.comboBox(
            units_box,
            self,
            "wavelength_unit",
            label="Wavelength:",
            items=WAVELENGTH_UNITS,
            sendSelectedValue=True,
            orientation=Qt.Horizontal,
            callback=self._units_changed,
        )
        gui.comboBox(
            units_box,
            self,
            "cd_unit",
            label="CD signal:",
            items=CD_SIGNAL_UNITS,
            sendSelectedValue=True,
            orientation=Qt.Horizontal,
            callback=self._units_changed,
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
            row,
            self,
            setting,
            label=label,
            orientation=Qt.Horizontal,
            callback=self.commit,
        )
        editor.setMinimumWidth(380)
        editor.setToolTip(str(getattr(self, setting)))
        editor.textChanged.connect(editor.setToolTip)
        gui.button(
            row,
            self,
            "Browse…",
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

    def _move_selected(self, step: int) -> None:
        """Move the selected files one place up (-1) or down (+1).

        The order of the list is the order of the titration points, so it
        must match the rows of the Titration Table.
        """
        items = [
            (path, index in self.selected_data_files)
            for index, path in enumerate(self.data_files)
        ]
        # Swap each selected item with an unselected neighbour, starting from
        # the end it is moving towards, so that adjacent selections move together.
        indices = range(1, len(items)) if step < 0 else range(len(items) - 2, -1, -1)
        changed = False
        for index in indices:
            neighbour = index + step
            if items[index][1] and not items[neighbour][1]:
                items[index], items[neighbour] = items[neighbour], items[index]
                changed = True
        if not changed:
            return
        self.data_files = [path for path, _ in items]
        self.data_file_names = self._file_labels()
        self.selected_data_files = [i for i, (_, chosen) in enumerate(items) if chosen]
        self._load_and_plot_raw_data()
        self.commit()

    def _clear_data_files(self) -> None:
        self.data_files = []
        self.cd_data = None
        self._refresh_file_list()
        self._refresh_plot()
        self.commit()

    @Inputs.titration
    def set_titration(self, table: Table | None) -> None:
        self.titration = table
        selected = list(self.selected_data_files)
        self.data_file_names = self._file_labels()
        self.selected_data_files = selected
        self.commit()

    def _points_available(self) -> bool:
        """Whether every data file can be given a row of the Titration Table."""
        return self.titration is not None and len(self.titration) == len(
            self.data_files
        )

    def _file_labels(self) -> list[str]:
        """List entries: the file names, numbered when a Titration Table is given."""
        names = [Path(path).name for path in self.data_files]
        if not self._points_available():
            return names
        try:
            ratios = self.titration.get_column("ratio")
        except (KeyError, ValueError):
            ratios = None
        labels = []
        for index, name in enumerate(names):
            label = f"{index + 1}. {name}"
            if ratios is not None:
                label += f"  (ratio {ratios[index]:g})"
            labels.append(label)
        return labels

    def _refresh_file_list(self) -> None:
        self.data_file_names = self._file_labels()
        # The selection only picks files to move; every file is always used.
        self.selected_data_files = []

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
                    x,
                    self.cd_data[name].to_numpy(dtype=float),
                    pen=pg.mkPen(pg.intColor(index, count), width=1.2),
                )
        self.plot_widget.enableAutoRange()

    def commit(self) -> None:
        """Send every loaded spectrum, uncorrected, as one table.

        A reference file becomes a ``raw_data`` row whose sample and role are
        the reference (``buffer``, ``sol_A`` or ``sol_B``), and a data file a
        ``raw_data`` row of role ``titration`` named after the file.
        """
        self.Error.load_failed.clear()
        self.Warning.cropped_wavelengths.clear()
        self.Warning.point_count_mismatch.clear()
        references = {
            ROLE_SOL_A: self.solution_a_file,
            ROLE_SOL_B: self.solution_b_file,
            ROLE_BUFFER: self.buffer_file,
        }
        try:
            series: dict[tuple[str, str], pd.Series] = {}
            for role, path in references.items():
                if path:
                    series[(role, STAGE_RAW, role)] = average_cd_series(path)
            for path in self.data_files:
                key = (Path(path).name, STAGE_RAW, ROLE_TITRATION)
                if key in series:
                    raise ValueError(f"Duplicate file name: {key[0]}")
                series[key] = average_cd_series(path)
            if not series:
                self.Outputs.data.send(None)
                return
            frame = pd.DataFrame(series)
            if frame.isna().any().any():
                # Keep only the wavelengths present in every file.
                frame = frame.dropna()
                if frame.empty:
                    raise ValueError("The selected files share no wavelengths")
                self.Warning.cropped_wavelengths(
                    count=len(frame),
                    low=float(frame.index.min()),
                    high=float(frame.index.max()),
                )
            points = None
            if self.titration is not None:
                if self._points_available():
                    # Number each data file by its place in the list, which is
                    # its row in the Titration Table. Later widgets pair the two
                    # by this number rather than by the order of the rows.
                    numbers = {
                        (Path(path).name, STAGE_RAW, ROLE_TITRATION): index + 1
                        for index, path in enumerate(self.data_files)
                    }
                    points = [numbers.get(column, np.nan) for column in frame.columns]
                else:
                    self.Warning.point_count_mismatch(
                        files=len(self.data_files), rows=len(self.titration)
                    )
            table = dataframe_to_orange_table(
                frame, self.wavelength_unit, self.cd_unit, points
            )
        except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
            self.Error.load_failed(str(exc))
            self.Outputs.data.send(None)
            return
        self.Outputs.data.send(table)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWCDDataLoader).run()
