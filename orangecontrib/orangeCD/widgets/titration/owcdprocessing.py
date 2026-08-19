"""Read and process circular-dichroism titration CSV files in Orange."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyqtgraph as pg
from Orange.data import ContinuousVariable, Domain, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget


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


def orange_table_to_dataframe(table: Table) -> pd.DataFrame:
    variables = (
        *table.domain.attributes,
        *table.domain.class_vars,
        *table.domain.metas,
    )
    return pd.DataFrame(
        {variable.name: table.get_column(variable) for variable in variables}
    )


def dataframe_to_orange_table(dataframe: pd.DataFrame) -> Table:
    output = dataframe.copy()
    output.columns = [
        " | ".join(map(str, col)) if isinstance(col, tuple) else str(col)
        for col in output.columns
    ]
    output.index.name = "Wavelength"
    output = output.reset_index()
    attributes = [col for col in output.columns if col != "Wavelength"]
    domain = Domain(
        [ContinuousVariable(col) for col in attributes],
        metas=[ContinuousVariable("Wavelength")],
    )
    return Table.from_numpy(
        domain,
        output[attributes].to_numpy(dtype=float),
        metas=output[["Wavelength"]].to_numpy(dtype=float),
    )


class OWCDTitrationProcessing(OWWidget):
    name = "CD Titration Processing"
    description = "Read CD spectra and apply titration reference corrections."
    icon = "icons/Titration.svg"
    priority = 20
    want_main_area = True
    resizing_enabled = True

    solution_a_file = Setting("")
    solution_b_file = Setting("")
    buffer_file = Setting("")
    data_files = Setting([])
    selected_data_files = Setting([])
    zero_range_min = Setting(440.0)
    zero_range_max = Setting(450.0)
    data_file_names: list[str] = []

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
        self.cd_data: pd.DataFrame | None = None
        self._updating_lines = False
        self._build_controls()
        self._build_plot()
        self._refresh_file_list()
        if self.data_files:
            self._load_and_plot_raw_data()

    def _build_controls(self) -> None:
        references = gui.widgetBox(self.controlArea, "Reference files")
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

        zero_box = gui.widgetBox(self.controlArea, "Zero-level wavelength range")
        gui.doubleSpin(
            zero_box, self, "zero_range_min", -1e6, 1e6,
            step=1.0, decimals=2, label="Lower limit",
            orientation="horizontal", callback=self._range_controls_changed,
        )
        gui.doubleSpin(
            zero_box, self, "zero_range_max", -1e6, 1e6,
            step=1.0, decimals=2, label="Upper limit",
            orientation="horizontal", callback=self._range_controls_changed,
        )
        gui.widgetLabel(
            zero_box,
            "Drag the two vertical lines on the raw-spectra plot to adjust the range.",
        ).setWordWrap(True)

        gui.button(self.controlArea, self, "Process data", callback=self.process)
        gui.rubber(self.controlArea)

    def _add_reference_selector(self, parent, label: str, setting: str) -> None:
        row = gui.hBox(parent)
        editor = gui.lineEdit(
            row, self, setting, label=label, orientation="horizontal"
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
        gui.widgetLabel(plot_box, "Raw spectra and zero-level averaging range")
        self.plot_widget = pg.PlotWidget(plot_box)
        self.plot_widget.setLabel("bottom", "Wavelength")
        self.plot_widget.setLabel("left", "Circular dichroism")
        self.plot_widget.showGrid(x=True, y=True, alpha=0.2)
        plot_box.layout().addWidget(self.plot_widget)

        line_pen = pg.mkPen("#d62728", width=2)
        hover_pen = pg.mkPen("#ff7f0e", width=3)
        self.lower_line = pg.InfiniteLine(
            pos=self.zero_range_min, angle=90, movable=True,
            pen=line_pen, hoverPen=hover_pen, label="Lower {value:.2f}",
        )
        self.upper_line = pg.InfiniteLine(
            pos=self.zero_range_max, angle=90, movable=True,
            pen=line_pen, hoverPen=hover_pen, label="Upper {value:.2f}",
        )
        self.lower_line.sigPositionChangeFinished.connect(self._range_lines_changed)
        self.upper_line.sigPositionChangeFinished.connect(self._range_lines_changed)
        self.plot_widget.addItem(self.lower_line)
        self.plot_widget.addItem(self.upper_line)

    def _choose_reference(self, setting: str) -> None:
        filename, _ = gui.QtWidgets.QFileDialog.getOpenFileName(
            self, "Select reference CSV", "", "CSV files (*.csv);;All files (*)"
        )
        if filename:
            setattr(self, setting, filename)

    def _choose_data_files(self) -> None:
        filenames, _ = gui.QtWidgets.QFileDialog.getOpenFileNames(
            self, "Select titration CSV files", "", "CSV files (*.csv);;All files (*)"
        )
        if filenames:
            self.data_files = sorted(dict.fromkeys(filenames), key=natural_key)
            self._refresh_file_list()
            self._load_and_plot_raw_data()

    def _clear_data_files(self) -> None:
        self.data_files = []
        self.cd_data = None
        self._refresh_file_list()
        self._refresh_plot()
        self.Outputs.data.send(None)

    def _refresh_file_list(self) -> None:
        self.data_file_names = [Path(path).name for path in self.data_files]
        self.selected_data_files = list(range(len(self.data_file_names)))

    def _load_and_plot_raw_data(self) -> None:
        try:
            self.cd_data = pd.concat(
                [average_cd_series(path) for path in self.data_files], axis=1
            )
        except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
            self.cd_data = None
            self.Error.processing_failed(str(exc))
        self._ensure_range_overlaps_data()
        self._refresh_plot()

    def _ensure_range_overlaps_data(self) -> None:
        if self.cd_data is None or self.cd_data.empty:
            return
        low = float(np.nanmin(self.cd_data.index.to_numpy(dtype=float)))
        high = float(np.nanmax(self.cd_data.index.to_numpy(dtype=float)))
        if self.zero_range_max < low or self.zero_range_min > high:
            span = high - low
            self.zero_range_min = high - 0.1 * span
            self.zero_range_max = high
        self._set_line_positions()

    def _refresh_plot(self) -> None:
        self.plot_widget.clear()
        if self.cd_data is not None:
            colours = pg.intColor
            x = self.cd_data.index.to_numpy(dtype=float)
            for index, name in enumerate(self.cd_data.columns):
                self.plot_widget.plot(
                    x, self.cd_data[name].to_numpy(dtype=float),
                    pen=pg.mkPen(colours(index, max(len(self.cd_data.columns), 1)), width=1.2),
                )
        self.plot_widget.addItem(self.lower_line)
        self.plot_widget.addItem(self.upper_line)
        self._set_line_positions()
        self.plot_widget.enableAutoRange()

    def _set_line_positions(self) -> None:
        if not hasattr(self, "lower_line"):
            return
        self._updating_lines = True
        self.lower_line.setValue(self.zero_range_min)
        self.upper_line.setValue(self.zero_range_max)
        self._updating_lines = False

    def _range_controls_changed(self) -> None:
        self._set_line_positions()
        if self.titration is not None and self.data_files:
            self.process()

    def _range_lines_changed(self) -> None:
        if self._updating_lines:
            return
        self.zero_range_min = float(min(self.lower_line.value(), self.upper_line.value()))
        self.zero_range_max = float(max(self.lower_line.value(), self.upper_line.value()))
        self._set_line_positions()
        if self.titration is not None and self.data_files:
            self.process()

    def _zero_region_mean(self, series: pd.Series) -> float:
        low = min(self.zero_range_min, self.zero_range_max)
        high = max(self.zero_range_min, self.zero_range_max)
        wavelength = series.index.to_numpy(dtype=float)
        mask = (wavelength >= low) & (wavelength <= high)
        if not np.any(mask):
            raise ValueError(
                f"No wavelength points fall inside the zero-level range {low:g} to {high:g}"
            )
        return float(series.iloc[np.flatnonzero(mask)].mean())

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
            ref_zero = self._zero_region_mean(
                reference_data["sol_A_buffer_subtracted"]
            )
            reference_data["sol_A_buffer_subtracted_zeroed"] = (
                reference_data["sol_A_buffer_subtracted"] - ref_zero
            )

            if self.cd_data is None:
                self._load_and_plot_raw_data()
            if self.cd_data is None:
                raise ValueError("The experimental CSV files could not be loaded")
            cd_data = self.cd_data
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
                    buffer_subtraction - reference_data["sol_A_buffer_subtracted"]
                )
                subtract_fraction_sol_b = (
                    sol_a_subtraction
                    - reference_data["sol_B_buffer_subtracted"]
                    * float(titration.iloc[index]["normalised_molar_ratio"])
                )
                zero_level = (
                    subtract_fraction_sol_b - self._zero_region_mean(subtract_fraction_sol_b)
                )

                plus_sol_a = (
                    zero_level + reference_data["sol_A_buffer_subtracted_zeroed"]
                )
                processed.update({
                    (series_name, "raw_data"): raw.round(3),
                    (series_name, "buffer_subtraction"): buffer_subtraction.round(3),
                    (series_name, "sol_A_subtraction"): sol_a_subtraction.round(3),
                    (series_name, "subtract_frac_sol_B"): subtract_fraction_sol_b.round(3),
                    (series_name, "set_to_zero"): round(zero_level, 3),
                    (series_name, "plus_sol_A"): plus_sol_a.round(3),
                })

            output = pd.DataFrame(processed, index=cd_data.index)
            orange_output = dataframe_to_orange_table(output)
        except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
            self.Error.processing_failed(str(exc))
            self.Outputs.data.send(None)
            return
        self.Outputs.data.send(orange_output)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWCDTitrationProcessing).run()
