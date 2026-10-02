#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert selected corrected CD spectra to mean-residue delta epsilon."""

from __future__ import annotations

import numpy as np
from pint import Quantity

from Orange.data import Domain, StringVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from . import Q_
from .utils import (
    CONCENTRATION_KEY,
    DELTA_EPSILON_UNIT,
    MDEG_PER_DELTA_A,
    PATHLENGTH_KEY,
    SPECTRUM_UNIT_KEY,
    UNIT_META,
    WAVELENGTH_META,
    InvalidWavelength,
    SpectraError,
    matching_spectra,
    quantity_string,
    reference_spectrum,
    spectrum_units,
    split_series_name,
    stages as spectra_stages,
    unit_string,
    wavelengths as spectra_wavelengths,
)

PATHLENGTH_UNIT = "centimeter"
MOLECULAR_WEIGHT_UNIT = "gram / mole"
CONCENTRATION_COLUMN = "working_concentration_a"
DEFAULT_CORRECTED_SERIES = "plus_sol_A"
DEFAULT_SOLUTION_A_SERIES = "sol_A_buffer_subtracted_zeroed"
DELTA_EPSILON_SUFFIX = "delta_epsilon"
DEFAULT_PATHLENGTH_CM = 1.0
DEFAULT_MEAN_RESIDUE_MW = 113.0
DEFAULT_SOLUTION_A_MW = 1.0


def calculate_delta_epsilon(
    cd: Quantity,
    concentration: Quantity,
    pathlength: Quantity,
    mean_residue_molecular_weight: Quantity,
    solution_a_molecular_weight: Quantity,
) -> Quantity:
    """Convert corrected CD to mean-residue delta epsilon.

    Every argument carries its own unit; the result is in
    ``DELTA_EPSILON_UNIT``. Raises ``ValueError`` for non-positive
    parameters and ``pint.DimensionalityError`` (a ``TypeError``) for
    incompatible units.
    """
    parameters = {
        "Concentration": concentration,
        "Pathlength": pathlength,
        "Mean residue molecular weight": mean_residue_molecular_weight,
        "Solution A molecular weight": solution_a_molecular_weight,
    }
    for name, quantity in parameters.items():
        value = float(quantity.magnitude)
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be a finite value greater than zero")

    delta_a = (cd.to(MDEG_PER_DELTA_A.units) / MDEG_PER_DELTA_A).to("dimensionless")
    residue_ratio = (
        mean_residue_molecular_weight / solution_a_molecular_weight
    ).to("dimensionless")
    return (
        delta_a * residue_ratio / (concentration * pathlength)
    ).to(DELTA_EPSILON_UNIT)


class OWDeltaEpsilon(OWWidget):
    name = "Delta Epsilon"
    description = (
        "Retain selected corrected CD spectra and add converted mean-residue "
        "delta epsilon spectra."
    )
    icon = "icons/DeltaEpsilon.svg"
    priority = 30
    want_main_area = False
    resizing_enabled = False

    class Inputs:
        data = Input("Processed CD Data", Table)
        titration = Input("Titration Table", Table)

    class Outputs:
        data = Output("CD and Delta Epsilon Spectra", Table)

    pathlength_cm = Setting(DEFAULT_PATHLENGTH_CM)
    mean_residue_molecular_weight = Setting(DEFAULT_MEAN_RESIDUE_MW)
    solution_a_molecular_weight = Setting(DEFAULT_SOLUTION_A_MW)
    corrected_series = Setting(DEFAULT_CORRECTED_SERIES)
    solution_a_series = Setting(DEFAULT_SOLUTION_A_SERIES)
    auto_commit = Setting(True)

    class Error(OWWidget.Error):
        missing_wavelength = Msg(
            "Input does not contain a Wavelength meta attribute."
        )
        invalid_wavelength = Msg(
            "Input attribute names are not all valid wavelengths."
        )
        no_series_names = Msg(
            "Input features do not use the expected 'sample | series' names."
        )
        no_corrected_series = Msg(
            "No features match the selected corrected series '{}'."
        )
        no_solution_a_series = Msg(
            "No feature matches the selected Solution A series '{}'."
        )
        invalid_parameter = Msg("{}")

    class Warning(OWWidget.Warning):
        no_titration = Msg(
            "Connect a Titration Table to provide the Solution A concentration."
        )
        missing_concentration = Msg(
            f"The Titration Table input does not contain a '{CONCENTRATION_COLUMN}' column."
        )

    def __init__(self) -> None:
        super().__init__()
        self.pathlength_cm = self._positive_float(
            self.pathlength_cm, DEFAULT_PATHLENGTH_CM
        )
        self.mean_residue_molecular_weight = self._positive_float(
            self.mean_residue_molecular_weight, DEFAULT_MEAN_RESIDUE_MW
        )
        self.solution_a_molecular_weight = self._positive_float(
            self.solution_a_molecular_weight, DEFAULT_SOLUTION_A_MW
        )
        if not isinstance(self.corrected_series, str):
            self.corrected_series = DEFAULT_CORRECTED_SERIES
        if not isinstance(self.solution_a_series, str):
            self.solution_a_series = DEFAULT_SOLUTION_A_SERIES

        self.data: Table | None = None
        self.titration_data: Table | None = None
        self.available_series: list[str] = []
        self._updating_series = False
        self._build_controls()

    @staticmethod
    def _positive_float(value: object, default: float) -> float:
        try:
            value = float(value)
        except (TypeError, ValueError):
            return default
        return value if np.isfinite(value) and value > 0 else default

    def _build_controls(self) -> None:
        series_box = gui.widgetBox(self.controlArea, "Data series")
        self.corrected_combo = gui.comboBox(
            series_box, self, "corrected_series",
            label="Corrected titration series", items=[],
            sendSelectedValue=True, valueType=str,
            orientation="horizontal", callback=self._series_changed,
        )
        self.solution_a_combo = gui.comboBox(
            series_box, self, "solution_a_series",
            label="Zeroed Solution A series", items=[],
            sendSelectedValue=True, valueType=str,
            orientation="horizontal", callback=self._series_changed,
        )
        note = gui.widgetLabel(
            series_box,
            "The output includes the original selected CD(mdeg) features and "
            "new matching features suffixed with '_delta_epsilon'.",
        )
        note.setWordWrap(True)

        conversion_box = gui.widgetBox(self.controlArea, "Conversion parameters")
        self.concentration_label = gui.widgetLabel(
            conversion_box, "Solution A concentration: connect a Titration Table"
        )
        self.concentration_label.setWordWrap(True)
        gui.doubleSpin(
            conversion_box, self, "pathlength_cm", 1e-12, 1e6,
            step=0.1, decimals=6,
            label=f"Pathlength ({self._unit_symbol(PATHLENGTH_UNIT)})",
            orientation="horizontal", callback=self.commit.deferred,
        )
        gui.doubleSpin(
            conversion_box, self, "mean_residue_molecular_weight", 1e-12, 1e6,
            step=1.0, decimals=3,
            label=f"Mean residue molecular weight ({self._unit_symbol(MOLECULAR_WEIGHT_UNIT)})",
            orientation="horizontal", callback=self.commit.deferred,
        )
        gui.doubleSpin(
            conversion_box, self, "solution_a_molecular_weight", 1e-12, 1e12,
            step=100.0, decimals=3,
            label=f"Solution A molecular weight ({self._unit_symbol(MOLECULAR_WEIGHT_UNIT)})",
            orientation="horizontal", callback=self.commit.deferred,
        )
        equation = gui.widgetLabel(
            conversion_box,
            "Delta epsilon = (CD / 32980 mdeg) x mean residue molecular weight / "
            "(concentration x pathlength x Solution A MW). "
            "Solution A concentration is taken from the connected Titration Table.",
        )
        equation.setWordWrap(True)
        gui.auto_commit(
            self.buttonsArea, self, "auto_commit", "Apply", commit=self.commit
        )

    @staticmethod
    def _unit_symbol(unit_name: str) -> str:
        return f"{Q_(1, unit_name).units:~}"

    @Inputs.data
    def set_data(self, data: Table | None) -> None:
        self.data = data
        self.Error.clear()
        self._update_series_controls()
        self.commit.now()

    @Inputs.titration
    def set_titration(self, table: Table | None) -> None:
        self.titration_data = table
        self.commit.deferred()

    def _solution_a_concentration(self) -> Quantity | None:
        if self.titration_data is None:
            return None
        try:
            variable = self.titration_data.domain[CONCENTRATION_COLUMN]
        except KeyError:
            return None
        unit = variable.attributes.get("unit")
        if not unit:
            return None
        values = self.titration_data.get_column(variable)
        if len(values) == 0 or not np.isfinite(values[0]):
            return None
        return Q_(float(values[0]), unit)

    def _update_series_controls(self) -> None:
        stages: list[str] = []
        if self.data is not None:
            try:
                stages = spectra_stages(self.data)
            except SpectraError:
                pass
        self.available_series = stages
        corrected = self._preferred(stages, self.corrected_series, DEFAULT_CORRECTED_SERIES)
        solution_a = self._preferred(stages, self.solution_a_series, DEFAULT_SOLUTION_A_SERIES)
        self._updating_series = True
        for combo, selected in (
            (self.corrected_combo, corrected),
            (self.solution_a_combo, solution_a),
        ):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(stages)
            if selected:
                combo.setCurrentText(selected)
            combo.blockSignals(False)
        self.corrected_series, self.solution_a_series = corrected, solution_a
        self._updating_series = False

    @staticmethod
    def _preferred(values: list[str], current: str, default: str) -> str:
        if default in values:
            return default
        if current in values:
            return current
        return values[0] if values else ""

    def _series_changed(self) -> None:
        if not self._updating_series:
            self.commit.deferred()

    def _wavelength_domain(self) -> Domain | None:
        """Return the input domain if it has a valid wavelength axis."""
        if self.data is None:
            return None
        try:
            spectra_wavelengths(self.data)
        except InvalidWavelength:
            self.Error.invalid_wavelength()
            return None
        except SpectraError:
            self.Error.missing_wavelength()
            return None
        return self.data.domain

    def _variables_matching(self, stage: str) -> list[tuple[int, str]]:
        if self.data is None or not stage:
            return []
        try:
            return matching_spectra(self.data, stage)
        except SpectraError:
            return []

    def _solution_a_variable(self) -> tuple[int, str] | None:
        return reference_spectrum(
            self._variables_matching(self.solution_a_series)
        )

    @staticmethod
    def _delta_name(variable: str) -> str:
        parsed = split_series_name(variable)
        if parsed is None:
            return f"{variable}_{DELTA_EPSILON_SUFFIX}"
        sample, stage = parsed
        return f"{sample} | {stage}_{DELTA_EPSILON_SUFFIX}"

    @gui.deferred
    def commit(self) -> None:
        self.Error.clear()
        self.Warning.clear()

        concentration = self._solution_a_concentration()
        if self.titration_data is None:
            self.Warning.no_titration()
            self.concentration_label.setText(
                "Solution A concentration: connect a Titration Table"
            )
        elif concentration is None:
            self.Warning.missing_concentration()
            self.concentration_label.setText(
                "Solution A concentration: unavailable from the connected Titration Table"
            )
        else:
            self.concentration_label.setText(
                f"Solution A concentration: {concentration:g} "
                f"(from Titration Table)"
            )

        if self.data is None:
            self.Outputs.data.send(None)
            return

        wavelength = self._wavelength_domain()
        if wavelength is None:
            self.Outputs.data.send(None)
            return
        if not self.available_series:
            self.Error.no_series_names()
            self.Outputs.data.send(None)
            return

        indexed_corrected = self._variables_matching(self.corrected_series)
        if not indexed_corrected:
            self.Error.no_corrected_series(self.corrected_series)
            self.Outputs.data.send(None)
            return
        solution_a = self._solution_a_variable()
        if solution_a is None:
            self.Error.no_solution_a_series(self.solution_a_series)
            self.Outputs.data.send(None)
            return

        if concentration is None:
            self.Outputs.data.send(None)
            return

        row_units = spectrum_units(self.data)
        raw_rows = [solution_a, *indexed_corrected]
        raw_variables_idx = [row for row, _ in raw_rows]
        raw_variable_names = [name for _, name in raw_rows]
        raw_units = [row_units[row] for row in raw_variables_idx]
        if not all(raw_units):
            self.Error.invalid_parameter(
                f"The input table has no '{SPECTRUM_UNIT_KEY}' attribute "
                f"or '{UNIT_META}' meta for the selected spectra."
            )
            self.Outputs.data.send(None)
            return

        pathlength = Q_(self.pathlength_cm, PATHLENGTH_UNIT)
        mean_residue_mw = Q_(
            self.mean_residue_molecular_weight, MOLECULAR_WEIGHT_UNIT
        )
        solution_a_mw = Q_(self.solution_a_molecular_weight, MOLECULAR_WEIGHT_UNIT)
        try:
            # (n_spectra, n_wavelengths), each row in its own unit.
            raw_values = self.data.X[raw_variables_idx]
            converted = np.vstack([
                calculate_delta_epsilon(
                    Q_(row_values, unit),
                    concentration,
                    pathlength,
                    mean_residue_mw,
                    solution_a_mw,
                ).magnitude
                for row_values, unit in zip(raw_values, raw_units)
            ])
        except (TypeError, ValueError) as exc:
            self.Error.invalid_parameter(str(exc))
            self.Outputs.data.send(None)
            return

        delta_unit = unit_string(DELTA_EPSILON_UNIT)
        delta_names = [self._delta_name(name) for name in raw_variable_names]
        output_domain = Domain(
            wavelength.attributes,
            metas=[
                next(m for m in wavelength.metas if m.name == WAVELENGTH_META),
                StringVariable(UNIT_META),
            ],
        )
        output = Table.from_numpy(
            output_domain,
            np.vstack((raw_values, converted)),
            metas=np.column_stack((
                np.asarray([*raw_variable_names, *delta_names], dtype=object),
                np.asarray(
                    [*raw_units, *[delta_unit] * len(delta_names)], dtype=object
                ),
            )),
            attributes=dict(self.data.attributes),
        )
        output.name = (
            f"{self.data.name} - CD and delta epsilon"
            if self.data.name else "CD and delta epsilon"
        )
        output.attributes.update({
            "corrected_series": self.corrected_series,
            "solution_a_series": self.solution_a_series,
            "delta_epsilon_suffix": DELTA_EPSILON_SUFFIX,
            CONCENTRATION_KEY: quantity_string(concentration),
            PATHLENGTH_KEY: quantity_string(pathlength),
            "mean_residue_molecular_weight": quantity_string(mean_residue_mw),
            "solution_a_molecular_weight": quantity_string(solution_a_mw),
        })
        self.Outputs.data.send(output)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWDeltaEpsilon).run()
