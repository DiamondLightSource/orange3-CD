#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert selected corrected CD spectra to mean-residue delta epsilon."""

from __future__ import annotations

import numpy as np

from Orange.data import ContinuousVariable, Domain, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget

MDEG_PER_DELTA_A = 32980.0
CONCENTRATION_UNITS = ("uM", "mM", "M")
CONCENTRATION_FACTORS = {"uM": 1e-6, "mM": 1e-3, "M": 1.0}
DEFAULT_CORRECTED_SERIES = "plus_sol_A"
DEFAULT_SOLUTION_A_SERIES = "sol_A_buffer_subtracted_zeroed"
DELTA_EPSILON_SUFFIX = "delta_epsilon"
DEFAULT_CONCENTRATION_UM = 19.659
DEFAULT_PATHLENGTH_CM = 1.0
DEFAULT_MEAN_RESIDUE_MW = 113.0
DEFAULT_SOLUTION_A_MW = 1.0


def split_series_name(variable_name: str) -> tuple[str, str] | None:
    if " | " not in variable_name:
        return None
    sample, stage = variable_name.rsplit(" | ", maxsplit=1)
    sample, stage = sample.strip(), stage.strip()
    return (sample, stage) if sample and stage else None


def calculate_delta_epsilon(
    cd_mdeg: np.ndarray,
    concentration: float,
    pathlength_cm: float,
    mean_residue_molecular_weight: float,
    solution_a_molecular_weight: float,
    concentration_unit: str,
) -> np.ndarray:
    """Convert corrected CD in mdeg to mean-residue delta epsilon."""
    if concentration_unit not in CONCENTRATION_FACTORS:
        raise ValueError(
            "Concentration unit must be one of: "
            + ", ".join(CONCENTRATION_UNITS)
        )
    parameters = {
        "Concentration": concentration,
        "Pathlength": pathlength_cm,
        "Mean residue molecular weight": mean_residue_molecular_weight,
        "Solution A molecular weight": solution_a_molecular_weight,
    }
    for name, value in parameters.items():
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be a finite value greater than zero")

    concentration_m = concentration * CONCENTRATION_FACTORS[concentration_unit]
    return (
        np.asarray(cd_mdeg, dtype=float)
        * mean_residue_molecular_weight
        / (
            MDEG_PER_DELTA_A
            * concentration_m
            * pathlength_cm
            * solution_a_molecular_weight
        )
    )


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

    class Outputs:
        data = Output("CD and Delta Epsilon Spectra", Table)

    concentration = Setting(DEFAULT_CONCENTRATION_UM)
    concentration_unit = Setting("uM")
    pathlength_cm = Setting(DEFAULT_PATHLENGTH_CM)
    mean_residue_molecular_weight = Setting(DEFAULT_MEAN_RESIDUE_MW)
    solution_a_molecular_weight = Setting(DEFAULT_SOLUTION_A_MW)
    corrected_series = Setting(DEFAULT_CORRECTED_SERIES)
    solution_a_series = Setting(DEFAULT_SOLUTION_A_SERIES)
    auto_commit = Setting(True)

    class Error(OWWidget.Error):
        missing_wavelength = Msg(
            "Input does not contain a continuous Wavelength meta attribute."
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

    def __init__(self) -> None:
        super().__init__()
        self.concentration = self._positive_float(
            self.concentration, DEFAULT_CONCENTRATION_UM
        )
        self.pathlength_cm = self._positive_float(
            self.pathlength_cm, DEFAULT_PATHLENGTH_CM
        )
        self.mean_residue_molecular_weight = self._positive_float(
            self.mean_residue_molecular_weight, DEFAULT_MEAN_RESIDUE_MW
        )
        self.solution_a_molecular_weight = self._positive_float(
            self.solution_a_molecular_weight, DEFAULT_SOLUTION_A_MW
        )
        self.concentration_unit = self._normalise_unit(self.concentration_unit)
        if not isinstance(self.corrected_series, str):
            self.corrected_series = DEFAULT_CORRECTED_SERIES
        if not isinstance(self.solution_a_series, str):
            self.solution_a_series = DEFAULT_SOLUTION_A_SERIES

        self.data: Table | None = None
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

    @staticmethod
    def _normalise_unit(value: object) -> str:
        if isinstance(value, int) and not isinstance(value, bool):
            return CONCENTRATION_UNITS[value] if 0 <= value < 3 else "uM"
        return value if value in CONCENTRATION_UNITS else "uM"

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
        gui.doubleSpin(
            conversion_box, self, "concentration", 1e-12, 1e12,
            step=0.1, decimals=6, label="Solution A concentration",
            orientation="horizontal", callback=self.commit.deferred,
        )
        gui.comboBox(
            conversion_box, self, "concentration_unit",
            label="Concentration unit", items=CONCENTRATION_UNITS,
            sendSelectedValue=True, valueType=str,
            orientation="horizontal", callback=self.commit.deferred,
        )
        gui.doubleSpin(
            conversion_box, self, "pathlength_cm", 1e-12, 1e6,
            step=0.1, decimals=6, label="Pathlength (cm)",
            orientation="horizontal", callback=self.commit.deferred,
        )
        gui.doubleSpin(
            conversion_box, self, "mean_residue_molecular_weight", 1e-12, 1e6,
            step=1.0, decimals=3, label="Mean residue molecular weight",
            orientation="horizontal", callback=self.commit.deferred,
        )
        gui.doubleSpin(
            conversion_box, self, "solution_a_molecular_weight", 1e-12, 1e12,
            step=100.0, decimals=3, label="Solution A molecular weight",
            orientation="horizontal", callback=self.commit.deferred,
        )
        equation = gui.widgetLabel(
            conversion_box,
            "Delta epsilon = CD(mdeg) x mean residue molecular weight / "
            "[32980 x concentration(M) x pathlength(cm) x Solution A MW]",
        )
        equation.setWordWrap(True)
        gui.auto_commit(
            self.buttonsArea, self, "auto_commit", "Apply", commit=self.commit
        )

    @Inputs.data
    def set_data(self, data: Table | None) -> None:
        self.data = data
        self.Error.clear()
        self._update_series_controls()
        self.commit.now()

    def _update_series_controls(self) -> None:
        stages: list[str] = []
        if self.data is not None:
            for variable in self.data.domain.attributes:
                parsed = split_series_name(variable.name)
                if parsed and parsed[1] not in stages:
                    stages.append(parsed[1])
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

    def _wavelength_variable(self) -> ContinuousVariable | None:
        if self.data is None:
            return None
        try:
            variable = self.data.domain["Wavelength"]
        except KeyError:
            return None
        return variable if (
            variable in self.data.domain.metas
            and isinstance(variable, ContinuousVariable)
        ) else None

    def _variables_matching(self, stage: str) -> list[ContinuousVariable]:
        if self.data is None or not stage:
            return []
        suffix = f" | {stage}"
        return [
            variable for variable in self.data.domain.attributes
            if isinstance(variable, ContinuousVariable)
            and variable.name.endswith(suffix)
        ]

    def _solution_a_variable(self) -> ContinuousVariable | None:
        candidates = self._variables_matching(self.solution_a_series)
        background = [
            variable for variable in candidates
            if variable.name.startswith("Background | ")
        ]
        return background[0] if background else (candidates[0] if candidates else None)

    @staticmethod
    def _delta_name(variable: ContinuousVariable) -> str:
        parsed = split_series_name(variable.name)
        if parsed is None:
            return f"{variable.name}_{DELTA_EPSILON_SUFFIX}"
        sample, stage = parsed
        return f"{sample} | {stage}_{DELTA_EPSILON_SUFFIX}"

    @gui.deferred
    def commit(self) -> None:
        self.Error.clear()
        if self.data is None:
            self.Outputs.data.send(None)
            return

        wavelength = self._wavelength_variable()
        if wavelength is None:
            self.Error.missing_wavelength()
            self.Outputs.data.send(None)
            return
        if not self.available_series:
            self.Error.no_series_names()
            self.Outputs.data.send(None)
            return

        corrected = self._variables_matching(self.corrected_series)
        if not corrected:
            self.Error.no_corrected_series(self.corrected_series)
            self.Outputs.data.send(None)
            return
        solution_a = self._solution_a_variable()
        if solution_a is None:
            self.Error.no_solution_a_series(self.solution_a_series)
            self.Outputs.data.send(None)
            return

        raw_variables = [solution_a, *corrected]
        raw_values = np.column_stack(
            [self.data.get_column(variable) for variable in raw_variables]
        )
        try:
            converted = calculate_delta_epsilon(
                raw_values,
                self.concentration,
                self.pathlength_cm,
                self.mean_residue_molecular_weight,
                self.solution_a_molecular_weight,
                self.concentration_unit,
            )
        except (TypeError, ValueError) as exc:
            self.Error.invalid_parameter(str(exc))
            self.Outputs.data.send(None)
            return

        delta_variables = [
            ContinuousVariable(self._delta_name(variable))
            for variable in raw_variables
        ]
        output_domain = Domain(
            [*raw_variables, *delta_variables], metas=[wavelength]
        )
        output = Table.from_numpy(
            output_domain,
            np.column_stack((raw_values, converted)),
            metas=np.asarray(
                self.data.get_column(wavelength), dtype=float
            ).reshape(-1, 1),
            ids=self.data.ids,
            attributes=dict(self.data.attributes),
        )
        output.name = (
            f"{self.data.name} - CD and delta epsilon"
            if self.data.name else "CD and delta epsilon"
        )
        output.attributes.update({
            "y_units": "mdeg and M^-1 cm^-1",
            "corrected_series": self.corrected_series,
            "solution_a_series": self.solution_a_series,
            "delta_epsilon_suffix": DELTA_EPSILON_SUFFIX,
            "solution_a_concentration_M": (
                self.concentration * CONCENTRATION_FACTORS[self.concentration_unit]
            ),
            "pathlength_cm": self.pathlength_cm,
            "mean_residue_molecular_weight": self.mean_residue_molecular_weight,
            "solution_a_molecular_weight": self.solution_a_molecular_weight,
        })
        self.Outputs.data.send(output)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWDeltaEpsilon).run()
