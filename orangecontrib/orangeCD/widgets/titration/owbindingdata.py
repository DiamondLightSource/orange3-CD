#!/usr/bin/env python3
"""Construct the CD Apps Binding data and Origin data at one wavelength."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from AnyQt.QtCore import Qt
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from ...units import Q_
from .utils import (
    CONCENTRATION_KEY,
    DELTA_EPSILON_UNIT,
    MDEG_PER_DELTA_A,
    MEASUREMENT_WAVELENGTH_KEY,
    PATHLENGTH_KEY,
    WAVELENGTH_UNIT_KEY,
    InvalidWavelength,
    SpectraError,
    matching_spectra,
    quantity_string,
    reference_spectrum,
    shared_unit,
    spectrum_units,
    table_quantity,
    table_unit,
    unit_string,
    unit_symbol,
)
from .utils import (
    stages as spectra_stages,
)
from .utils import (
    wavelengths as spectra_wavelengths,
)

# Preferred series, most-corrected first; used when nothing valid is saved.
DEFAULT_DATA_SERIES = ("plus_sol_A_delta_epsilon", "raw_data")
DEFAULT_SOLUTION_A_SERIES = ("sol_A_buffer_subtracted_delta_epsilon", "sol_A")
# Origin/CD Apps convention: concentration [B] is reported in molar.
CONCENTRATION_UNIT = "molar"


def continuous_column(table: Table, names: tuple[str, ...]) -> np.ndarray | None:
    """Return the first available continuous column matching *names*."""
    for name in names:
        try:
            variable = table.domain[name]
        except KeyError:
            continue
        if isinstance(variable, ContinuousVariable):
            return np.asarray(table.get_column(variable), dtype=float)
    return None


class OWBindingData(OWWidget):
    name = "Binding Data"
    description = (
        "Construct the Binding data and Origin data columns at a selected wavelength."
    )
    icon = "icons/BindingData.svg"
    priority = 40
    want_main_area = True
    resizing_enabled = True

    class Inputs:
        spectra = Input("CD and Delta Epsilon Spectra", Table)
        titration = Input("Titration Table", Table)

    class Outputs:
        data = Output("Binding Data", Table)

    wavelength = Setting(400.0)
    data_series = Setting("")
    solution_a_series = Setting("")
    absolute_change = Setting(True)

    class Error(OWWidget.Error):
        missing_wavelength = Msg("Input has no 'Spectrum' string meta naming the rows.")
        invalid_wavelength = Msg("Input attribute names are not all valid wavelengths.")
        missing_solution_a = Msg("No Solution A CD feature matching '{}' was found.")
        missing_data = Msg("No titration data features matching '{}' were found.")
        missing_ratios = Msg("Titration input has no ratio column.")
        missing_conversion_metadata = Msg(
            "The spectra table is missing '{}' conversion metadata."
        )
        invalid_conversion_metadata = Msg("{} must be greater than zero.")
        mismatched_points = Msg("{}")
        incompatible_units = Msg("{}")

    def __init__(self) -> None:
        super().__init__()
        self.spectra: Table | None = None
        self.titration: Table | None = None
        self._wavelengths: np.ndarray | None = None
        self._wavelength_error = None
        self._moving_line = False
        self._updating_series = False
        self._build_controls()
        self._build_plot()

    def _wavelength_unit(self) -> str | None:
        if self.spectra is None:
            return None
        return table_unit(self.spectra, WAVELENGTH_UNIT_KEY)

    def _plotted_cd_unit(self) -> str | None:
        """The CD unit shared by the plotted spectra, if there is one."""
        if self.spectra is None:
            return None
        row_units = spectrum_units(self.spectra)
        rows = [row for row, _ in self._matching(self.data_series)]
        reference = self._reference()
        if reference is not None:
            rows.append(reference[0])
        return shared_unit(row_units[row] for row in rows) if rows else None

    def _build_controls(self) -> None:
        box = gui.widgetBox(self.controlArea, "Binding-data settings")
        self.wavelength_spin = gui.doubleSpin(
            box,
            self,
            "wavelength",
            -1e6,
            1e6,
            step=1.0,
            decimals=3,
            label="Measurement wavelength",
            orientation=Qt.Horizontal,
            callback=self._wavelength_control_changed,
        )
        self.data_combo = gui.comboBox(
            box,
            self,
            "data_series",
            label="Titration data series",
            items=[],
            sendSelectedValue=True,
            orientation=Qt.Horizontal,
            callback=self._series_changed,
        )
        self.solution_a_combo = gui.comboBox(
            box,
            self,
            "solution_a_series",
            label="Solution A series",
            items=[],
            sendSelectedValue=True,
            orientation=Qt.Horizontal,
            callback=self._series_changed,
        )
        gui.checkBox(
            box,
            self,
            "absolute_change",
            "Use absolute change from Solution A",
            callback=self.commit,
        )
        note = gui.widgetLabel(
            box,
            "Origin concentration [B] is calculated as titration point "
            "multiplied by the molar concentration of Solution A.",
        )
        note.setWordWrap(True)
        gui.rubber(self.controlArea)

    def _build_plot(self) -> None:
        box = gui.vBox(self.mainArea)
        gui.widgetLabel(box, "CD spectra and selected wavelength")
        self.plot = pg.PlotWidget(box)
        self.plot.setLabel("bottom", "Wavelength")
        self.plot.setLabel("left", "CD")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        box.layout().addWidget(self.plot)

        self.line = pg.InfiniteLine(
            pos=self.wavelength,
            angle=90,
            movable=True,
            pen=pg.mkPen("#d62728", width=2),
            hoverPen=pg.mkPen("#ff7f0e", width=3),
            label="{value:.3f}",
        )
        self.line.sigPositionChangeFinished.connect(self._line_changed)
        self.plot.addItem(self.line)

    @Inputs.spectra
    def set_spectra(self, data: Table | None) -> None:
        self.spectra = data
        self._update_series_controls()
        self._load_wavelengths()
        self._snap()
        self._refresh_plot()
        self.commit()

    def _update_series_controls(self) -> None:
        stages: list[str] = []
        if self.spectra is not None:
            try:
                stages = spectra_stages(self.spectra)
            except SpectraError:
                pass
        data_stage = self._preferred(stages, self.data_series, DEFAULT_DATA_SERIES)
        solution_a = self._preferred(
            stages, self.solution_a_series, DEFAULT_SOLUTION_A_SERIES
        )
        self._updating_series = True
        for combo, selected in (
            (self.data_combo, data_stage),
            (self.solution_a_combo, solution_a),
        ):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(stages)
            if selected:
                combo.setCurrentText(selected)
            combo.blockSignals(False)
        # Assigning to an empty combo makes Orange warn; keep the old value.
        if data_stage:
            self.data_series = data_stage
        if solution_a:
            self.solution_a_series = solution_a
        self._updating_series = False

    @staticmethod
    def _preferred(values: list[str], current: str, defaults: tuple[str, ...]) -> str:
        if current in values:
            return current
        for default in defaults:
            if default in values:
                return default
        return values[0] if values else ""

    @Inputs.titration
    def set_titration(self, data: Table | None) -> None:
        self.titration = data
        self.commit()

    def _load_wavelengths(self) -> None:
        self._wavelengths = None
        self._wavelength_error = None
        if self.spectra is None:
            return
        try:
            self._wavelengths = spectra_wavelengths(self.spectra)
        except InvalidWavelength:
            self._wavelength_error = self.Error.invalid_wavelength
        except SpectraError:
            self._wavelength_error = self.Error.missing_wavelength

    def _snap(self, requested: float | None = None) -> int | None:
        if self._wavelengths is None:
            return None
        target = self.wavelength if requested is None else requested
        index = int(np.argmin(np.abs(self._wavelengths - target)))
        self.wavelength = float(self._wavelengths[index])
        self._moving_line = True
        self.line.setValue(self.wavelength)
        self._moving_line = False
        return index

    def _wavelength_control_changed(self) -> None:
        self._snap()
        self.commit()

    def _line_changed(self) -> None:
        if not self._moving_line:
            self._snap(float(self.line.value()))
            self.commit()

    def _series_changed(self) -> None:
        if not self._updating_series:
            self._refresh_plot()
            self.commit()

    def _matching(self, stage: str) -> list[tuple[int, str]]:
        if self.spectra is None:
            return []
        try:
            return matching_spectra(self.spectra, stage)
        except SpectraError:
            return []

    def _reference(self) -> tuple[int, str] | None:
        return reference_spectrum(self._matching(self.solution_a_series))

    def _ratios(self) -> np.ndarray | None:
        if self.titration is None:
            return None
        return continuous_column(
            self.titration,
            ("ratio", "normalised_molar_ratio", "molar_ratio"),
        )

    def _metadata_quantity(self, name: str):
        if self.spectra is None or name not in self.spectra.attributes:
            self.Error.missing_conversion_metadata(name)
            return None
        try:
            return table_quantity(self.spectra, name)
        except ValueError:
            self.Error.invalid_conversion_metadata(name)
            return None

    def _residue_ratio(self) -> float:
        """Mean-residue / Solution A molecular weight used for delta epsilon."""
        try:
            mean_residue = table_quantity(self.spectra, "mean_residue_molecular_weight")
            solution_a = table_quantity(self.spectra, "solution_a_molecular_weight")
            return float((mean_residue / solution_a).to("dimensionless").magnitude)
        except (KeyError, ValueError, TypeError):
            return 1.0

    def _refresh_plot(self) -> None:
        self.plot.clear()
        wavelength_symbol = unit_symbol(self._wavelength_unit())
        self.plot.setLabel("bottom", "Wavelength", units=wavelength_symbol)
        self.plot.setLabel("left", "CD", units=unit_symbol(self._plotted_cd_unit()))
        self.line.label.setFormat(
            "{value:.3f}" + (f" {wavelength_symbol}" if wavelength_symbol else "")
        )
        self.wavelength_spin.setSuffix(
            f" {wavelength_symbol}" if wavelength_symbol else ""
        )
        if self.spectra is not None and self._wavelengths is not None:
            indexed_rows = self._matching(self.data_series)
            variable_indices = [i[0] for i in indexed_rows]
            reference = self._reference()
            if reference is not None:
                variable_indices.insert(0, reference[0])
            for index, variable in enumerate(variable_indices):
                self.plot.plot(
                    self._wavelengths,
                    np.asarray(self.spectra.X[variable], dtype=float),
                    pen=pg.mkPen(
                        pg.intColor(index, max(len(variable_indices), 1)),
                        width=1.2,
                    ),
                )
        self.plot.addItem(self.line)
        self._snap()
        self.plot.enableAutoRange()

    def commit(self) -> None:
        self.Error.clear()
        if self._wavelength_error is not None:
            self._wavelength_error()
            self.Outputs.data.send(None)
            return
        if self.spectra is None or self.titration is None:
            self.Outputs.data.send(None)
            return

        row = self._snap()
        if row is None:
            self.Outputs.data.send(None)
            return

        reference = self._reference()
        if reference is None:
            self.Error.missing_solution_a(self.solution_a_series)
            self.Outputs.data.send(None)
            return

        indexed_data = self._matching(self.data_series)
        if not indexed_data:
            self.Error.missing_data(self.data_series)
            self.Outputs.data.send(None)
            return

        ratios = self._ratios()
        if ratios is None:
            self.Error.missing_ratios()
            self.Outputs.data.send(None)
            return
        if len(ratios) != len(indexed_data):
            self.Error.mismatched_points(
                f"There are {len(ratios)} titration ratios but "
                f"{len(indexed_data)} spectra"
            )
            self.Outputs.data.send(None)
            return

        concentration = self._metadata_quantity(CONCENTRATION_KEY)
        pathlength = self._metadata_quantity(PATHLENGTH_KEY)
        if concentration is None or pathlength is None:
            self.Outputs.data.send(None)
            return

        variables_idx = [reference[0], *(i[0] for i in indexed_data)]
        sample_names = [reference[1], *(i[1] for i in indexed_data)]

        # Each spectrum keeps the unit it was given on input; bring them all
        # to the unit of the first row before taking differences.
        row_units = spectrum_units(self.spectra)
        units = [row_units[index] for index in variables_idx]
        if not all(units):
            self.Error.incompatible_units(
                "The spectra table has no spectrum unit "
                "(attribute 'spectrum_unit' or 'Unit' meta)."
            )
            self.Outputs.data.send(None)
            return
        cd_unit = units[0]
        try:
            cd = Q_(
                np.asarray(
                    [
                        Q_(float(self.spectra.X[index, row]), unit)
                        .to(cd_unit)
                        .magnitude
                        for index, unit in zip(variables_idx, units)
                    ]
                ),
                cd_unit,
            )
            # Either an angular CD signal or a delta epsilon spectrum can be
            # differenced; anything else cannot give a binding signal.
            is_delta_epsilon = cd.check(DELTA_EPSILON_UNIT)
            if not is_delta_epsilon:
                cd.to(MDEG_PER_DELTA_A.units)
        except TypeError as exc:
            self.Error.incompatible_units(
                "The selected series must be in CD units (e.g. millidegree) "
                f"or delta epsilon units: {exc}"
            )
            self.Outputs.data.send(None)
            return

        titration_point = np.concatenate(([0.0], np.asarray(ratios, dtype=float)))

        change = np.zeros_like(cd.magnitude) * cd.units
        if self.absolute_change:
            change[1:] = np.abs(cd[1:] - cd[0])
        else:
            change[1:] = cd[1:] - cd[0]

        if is_delta_epsilon:
            delta_epsilon = change.to(DELTA_EPSILON_UNIT)
            # Delta epsilon spectra are scaled to mean-residue values by the
            # Delta Epsilon widget; undo that to recover the measured delta A.
            delta_a = (
                delta_epsilon * concentration * pathlength / self._residue_ratio()
            ).to("dimensionless")
        else:
            delta_a = (change.to(MDEG_PER_DELTA_A.units) / MDEG_PER_DELTA_A).to(
                "dimensionless"
            )
            delta_epsilon = (delta_a / (concentration * pathlength)).to(
                DELTA_EPSILON_UNIT
            )
        binding_stoichiometry = titration_point / (titration_point + 1.0)

        # CD Apps Origin data column 1:
        # concentration [B] = titration point * host concentration.
        concentration_b = (titration_point * concentration).to(CONCENTRATION_UNIT)

        cd_variable = ContinuousVariable("CD")
        cd_variable.attributes["unit"] = unit_string(cd_unit)
        change_variable = ContinuousVariable("Change in CD")
        change_variable.attributes["unit"] = unit_string(cd_unit)
        delta_epsilon_variable = ContinuousVariable("Delta Epsilon")
        delta_epsilon_variable.attributes["unit"] = unit_string(DELTA_EPSILON_UNIT)
        concentration_b_variable = ContinuousVariable("Conc [B]")
        concentration_b_variable.attributes["unit"] = unit_string(CONCENTRATION_UNIT)
        domain = Domain(
            [
                ContinuousVariable("Titration point"),
                cd_variable,
                change_variable,
                ContinuousVariable("Delta A"),
                delta_epsilon_variable,
                ContinuousVariable("Binding Stoichiometry"),
                concentration_b_variable,
            ],
            metas=[StringVariable("Sample")],
        )

        output = Table.from_numpy(
            domain,
            np.column_stack(
                (
                    titration_point,
                    cd.magnitude,
                    change.magnitude,
                    delta_a.magnitude,
                    delta_epsilon.magnitude,
                    binding_stoichiometry,
                    concentration_b.magnitude,
                )
            ),
            metas=np.asarray(sample_names, dtype=object).reshape(-1, 1),
        )
        wavelength_unit = self._wavelength_unit()
        measurement = Q_(self.wavelength, wavelength_unit or "dimensionless")
        output.name = f"Binding and Origin data at {self.wavelength:g}" + (
            f" {unit_symbol(wavelength_unit)}" if wavelength_unit else ""
        )
        output.attributes.update(
            {
                MEASUREMENT_WAVELENGTH_KEY: quantity_string(measurement),
                "data_series": self.data_series,
                "solution_a_series": self.solution_a_series,
                "absolute_change": self.absolute_change,
                CONCENTRATION_KEY: quantity_string(concentration),
                PATHLENGTH_KEY: quantity_string(pathlength),
                "delta_epsilon_calculation": (
                    f"Delta_A / ({CONCENTRATION_KEY} * {PATHLENGTH_KEY})"
                ),
                "origin_concentration_b_calculation": (
                    f"Titration_point * {CONCENTRATION_KEY}"
                ),
            }
        )
        self.Outputs.data.send(output)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWBindingData).run()
