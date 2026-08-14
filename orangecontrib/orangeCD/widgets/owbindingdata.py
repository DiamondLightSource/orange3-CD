#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construct the CD Apps Binding data and Origin data at one wavelength."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget

DEFAULT_CORRECTED_SERIES = "plus_sol_A"
DEFAULT_SOLUTION_A_SERIES = "sol_A_buffer_subtracted_zeroed"
MDEG_PER_DELTA_A = 32980.0


def split_series_name(name: str) -> tuple[str, str] | None:
    """Split a feature name of the form ``sample | processing_stage``."""
    if " | " not in name:
        return None
    sample, stage = name.rsplit(" | ", maxsplit=1)
    sample, stage = sample.strip(), stage.strip()
    return (sample, stage) if sample and stage else None


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
        "Construct the Binding data and Origin data columns at a selected "
        "wavelength."
    )
    icon = "icons/Titration.svg"
    priority = 40
    want_main_area = True
    resizing_enabled = True

    class Inputs:
        spectra = Input("CD and Delta Epsilon Spectra", Table)
        titration = Input("Titration Table", Table)

    class Outputs:
        data = Output("Binding Data", Table)

    wavelength = Setting(355.0)
    corrected_series = Setting(DEFAULT_CORRECTED_SERIES)
    solution_a_series = Setting(DEFAULT_SOLUTION_A_SERIES)
    absolute_change = Setting(True)

    class Error(OWWidget.Error):
        missing_wavelength = Msg("Input has no continuous Wavelength meta.")
        missing_solution_a = Msg(
            "No Solution A CD feature matching '{}' was found."
        )
        missing_corrected = Msg(
            "No corrected CD features matching '{}' were found."
        )
        missing_ratios = Msg("Titration input has no ratio column.")
        missing_conversion_metadata = Msg(
            "The spectra table is missing '{}' conversion metadata."
        )
        invalid_conversion_metadata = Msg("{} must be greater than zero.")
        mismatched_points = Msg("{}")

    def __init__(self) -> None:
        super().__init__()
        self.spectra: Table | None = None
        self.titration: Table | None = None
        self._wavelengths: np.ndarray | None = None
        self._moving_line = False
        self._build_controls()
        self._build_plot()

    def _build_controls(self) -> None:
        box = gui.widgetBox(self.controlArea, "Binding-data settings")
        gui.doubleSpin(
            box,
            self,
            "wavelength",
            -1e6,
            1e6,
            step=1.0,
            decimals=3,
            label="Measurement wavelength (nm)",
            orientation="horizontal",
            callback=self._wavelength_control_changed,
        )
        gui.lineEdit(
            box,
            self,
            "corrected_series",
            label="Corrected CD series",
            orientation="horizontal",
            callback=self._settings_changed,
        )
        gui.lineEdit(
            box,
            self,
            "solution_a_series",
            label="Zeroed Solution A series",
            orientation="horizontal",
            callback=self._settings_changed,
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
        gui.widgetLabel(box, "Corrected CD spectra and selected wavelength")
        self.plot = pg.PlotWidget(box)
        self.plot.setLabel("bottom", "Wavelength", units="nm")
        self.plot.setLabel("left", "CD", units="mdeg")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        box.layout().addWidget(self.plot)

        self.line = pg.InfiniteLine(
            pos=self.wavelength,
            angle=90,
            movable=True,
            pen=pg.mkPen("#d62728", width=2),
            hoverPen=pg.mkPen("#ff7f0e", width=3),
            label="{value:.3f} nm",
        )
        self.line.sigPositionChangeFinished.connect(self._line_changed)
        self.plot.addItem(self.line)

    @Inputs.spectra
    def set_spectra(self, data: Table | None) -> None:
        self.spectra = data
        self._load_wavelengths()
        self._snap()
        self._refresh_plot()
        self.commit()

    @Inputs.titration
    def set_titration(self, data: Table | None) -> None:
        self.titration = data
        self.commit()

    def _load_wavelengths(self) -> None:
        self._wavelengths = None
        if self.spectra is None:
            return
        try:
            variable = self.spectra.domain["Wavelength"]
        except KeyError:
            return
        if (
            variable in self.spectra.domain.metas
            and isinstance(variable, ContinuousVariable)
        ):
            values = np.asarray(self.spectra.get_column(variable), dtype=float)
            if values.size and np.all(np.isfinite(values)):
                self._wavelengths = values

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

    def _settings_changed(self) -> None:
        self._refresh_plot()
        self.commit()

    def _matching(self, stage: str) -> list[ContinuousVariable]:
        if self.spectra is None:
            return []
        suffix = f" | {stage}"
        return [
            variable
            for variable in self.spectra.domain.attributes
            if isinstance(variable, ContinuousVariable)
            and variable.name.endswith(suffix)
        ]

    def _reference(self) -> ContinuousVariable | None:
        candidates = self._matching(self.solution_a_series)
        background = [
            variable
            for variable in candidates
            if variable.name.startswith("Background | ")
        ]
        return background[0] if background else (
            candidates[0] if candidates else None
        )

    def _ratios(self) -> np.ndarray | None:
        if self.titration is None:
            return None
        return continuous_column(
            self.titration,
            ("ratio", "normalised_molar_ratio", "molar_ratio"),
        )

    def _metadata_float(self, name: str) -> float | None:
        if self.spectra is None or name not in self.spectra.attributes:
            self.Error.missing_conversion_metadata(name)
            return None
        try:
            value = float(self.spectra.attributes[name])
        except (TypeError, ValueError):
            self.Error.invalid_conversion_metadata(name)
            return None
        if not np.isfinite(value) or value <= 0:
            self.Error.invalid_conversion_metadata(name)
            return None
        return value

    def _refresh_plot(self) -> None:
        self.plot.clear()
        if self.spectra is not None and self._wavelengths is not None:
            variables = self._matching(self.corrected_series)
            reference = self._reference()
            if reference is not None:
                variables.insert(0, reference)
            for index, variable in enumerate(variables):
                self.plot.plot(
                    self._wavelengths,
                    np.asarray(
                        self.spectra.get_column(variable), dtype=float
                    ),
                    pen=pg.mkPen(
                        pg.intColor(index, max(len(variables), 1)),
                        width=1.2,
                    ),
                )
        self.plot.addItem(self.line)
        self._snap()
        self.plot.enableAutoRange()

    def commit(self) -> None:
        self.Error.clear()
        if self.spectra is None or self.titration is None:
            self.Outputs.data.send(None)
            return

        row = self._snap()
        if row is None:
            self.Error.missing_wavelength()
            self.Outputs.data.send(None)
            return

        reference = self._reference()
        if reference is None:
            self.Error.missing_solution_a(self.solution_a_series)
            self.Outputs.data.send(None)
            return

        corrected = self._matching(self.corrected_series)
        if not corrected:
            self.Error.missing_corrected(self.corrected_series)
            self.Outputs.data.send(None)
            return

        ratios = self._ratios()
        if ratios is None:
            self.Error.missing_ratios()
            self.Outputs.data.send(None)
            return
        if len(ratios) != len(corrected):
            self.Error.mismatched_points(
                f"There are {len(ratios)} titration ratios but "
                f"{len(corrected)} spectra"
            )
            self.Outputs.data.send(None)
            return

        concentration_m = self._metadata_float(
            "solution_a_concentration_M"
        )
        pathlength_cm = self._metadata_float("pathlength_cm")
        if concentration_m is None or pathlength_cm is None:
            self.Outputs.data.send(None)
            return

        raw_variables = [reference, *corrected]
        cd_mdeg = np.asarray(
            [
                self.spectra.get_column(variable)[row]
                for variable in raw_variables
            ],
            dtype=float,
        )
        titration_point = np.concatenate(
            ([0.0], np.asarray(ratios, dtype=float))
        )

        change_mdeg = np.zeros_like(cd_mdeg)
        if self.absolute_change:
            change_mdeg[1:] = np.abs(cd_mdeg[1:] - cd_mdeg[0])
        else:
            change_mdeg[1:] = cd_mdeg[1:] - cd_mdeg[0]

        delta_a = change_mdeg / MDEG_PER_DELTA_A
        delta_epsilon = delta_a / (concentration_m * pathlength_cm)
        binding_stoichiometry = titration_point / (titration_point + 1.0)

        # CD Apps Origin data column 1:
        # concentration [B] (M) = titration point * host concentration (M).
        concentration_b_m = titration_point * concentration_m

        sample_names = ["Solution A"] + [
            split_series_name(variable.name)[0] for variable in corrected
        ]
        domain = Domain(
            [
                ContinuousVariable("Titration point"),
                ContinuousVariable("CD(mdeg)"),
                ContinuousVariable("Change in CD(mdeg)"),
                ContinuousVariable("Delta A"),
                ContinuousVariable("Delta Epsilon"),
                ContinuousVariable("Binding Stoichiometry"),
                ContinuousVariable("Conc [B] (Molar)"),
            ],
            metas=[StringVariable("Sample")],
        )
        output = Table.from_numpy(
            domain,
            np.column_stack(
                (
                    titration_point,
                    cd_mdeg,
                    change_mdeg,
                    delta_a,
                    delta_epsilon,
                    binding_stoichiometry,
                    concentration_b_m,
                )
            ),
            metas=np.asarray(sample_names, dtype=object).reshape(-1, 1),
        )
        output.name = f"Binding and Origin data at {self.wavelength:g} nm"
        output.attributes.update(
            {
                "measurement_wavelength_nm": self.wavelength,
                "corrected_series": self.corrected_series,
                "solution_a_series": self.solution_a_series,
                "absolute_change": self.absolute_change,
                "solution_a_concentration_M": concentration_m,
                "pathlength_cm": pathlength_cm,
                "delta_epsilon_calculation": (
                    "Delta_A / "
                    "(solution_a_concentration_M * pathlength_cm)"
                ),
                "origin_concentration_b_calculation": (
                    "Titration_point * solution_a_concentration_M"
                ),
            }
        )
        self.Outputs.data.send(output)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWBindingData).run()
