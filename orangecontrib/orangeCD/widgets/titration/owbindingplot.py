#!/usr/bin/env python3
"""Scatter plot and manual Hill-equation fitting for Binding Data tables."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from AnyQt.QtCore import Qt
from Orange.data import ContinuousVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from ...units import Q_
from .binding_models import MODELS, BindingModel
from .utils import MEASUREMENT_WAVELENGTH_KEY, table_quantity

# The Hill fit of CD Apps used [B] in micromolar against delta A per cm of
# pathlength. Tables without those columns fall back to the plain ones.
DEFAULT_X = "Hill Conc [B]"
DEFAULT_Y = "Hill Delta A"
FALLBACK_X = "Titration point"
FALLBACK_Y = "Delta A"


class OWBindingPlot(OWWidget):
    name = "Binding Plot"
    description = "Plot Binding Data and fit a Hill or BiHill binding model."
    icon = "icons/BindingPlot.svg"
    priority = 50
    want_main_area = True
    resizing_enabled = True

    class Inputs:
        data = Input("Binding Data", Table)

    class Outputs:
        fit_results = Output("Hill Fit Results", Table)

    x_variable = Setting(DEFAULT_X)
    y_variable = Setting(DEFAULT_Y)
    point_size = Setting(9)
    line_width = Setting(2.0)
    # Index into ``MODELS``.
    fit_model = Setting(0)

    class Error(OWWidget.Error):
        no_continuous_data = Msg("Input contains no continuous variables.")
        invalid_selection = Msg("Select valid x and y variables.")
        fit_failed = Msg("{}")

    def __init__(self) -> None:
        super().__init__()
        self.data: Table | None = None
        self.variable_names: list[str] = []
        self._updating_controls = False
        self._fit_curve: dict[np.ndarray, np.ndarray, np.ndarray] | None = None
        self._build_controls()
        self._build_plot()

    def _build_controls(self) -> None:
        variables_box = gui.widgetBox(self.controlArea, "Plot variables")
        self.x_combo = gui.comboBox(
            variables_box,
            self,
            "x_variable",
            label="X variable",
            items=[],
            sendSelectedValue=True,
            orientation=Qt.Horizontal,
            callback=self._selection_changed,
        )
        self.y_combo = gui.comboBox(
            variables_box,
            self,
            "y_variable",
            label="Y variable",
            items=[],
            sendSelectedValue=True,
            orientation=Qt.Horizontal,
            callback=self._selection_changed,
        )
        self.x_combo.setMinimumWidth(220)
        self.y_combo.setMinimumWidth(220)

        appearance_box = gui.widgetBox(self.controlArea, "Appearance")
        gui.spin(
            appearance_box,
            self,
            "point_size",
            2,
            30,
            label="Point size",
            orientation=Qt.Horizontal,
            callback=self._redraw,
        )
        gui.doubleSpin(
            appearance_box,
            self,
            "line_width",
            0.1,
            10.0,
            step=0.1,
            decimals=1,
            label="Fit line width",
            orientation=Qt.Horizontal,
            callback=self._redraw,
        )

        fit_box = gui.widgetBox(self.controlArea, "Hill fit")
        if not isinstance(self.fit_model, int) or not (
            0 <= self.fit_model < len(MODELS)
        ):
            self.fit_model = 0
        self.equation = gui.widgetLabel(fit_box, self.binding_model.equation)
        self.equation.setWordWrap(True)

        gui.comboBox(
            fit_box,
            self,
            "fit_model",
            label="Select fit model",
            items=[model.name for model in MODELS],
            sendSelectedValue=False,
            orientation=Qt.Horizontal,
            callback=self._select_fit_model_changed,
        )
        gui.button(
            fit_box,
            self,
            "Fit model",
            callback=self.fit,
        )
        gui.button(
            fit_box,
            self,
            "Clear fit",
            callback=self.clear_fit,
        )
        gui.rubber(self.controlArea)

    @staticmethod
    def _unit_symbol(unit_name: str | None) -> str | None:
        if not unit_name:
            return None
        return f"{Q_(1, unit_name).units:~}"

    def _measurement_wavelength(self) -> str | None:
        """The wavelength the binding data was taken at, e.g. ``"400 nm"``."""
        if self.data is None or MEASUREMENT_WAVELENGTH_KEY not in self.data.attributes:
            return None
        try:
            quantity = table_quantity(self.data, MEASUREMENT_WAVELENGTH_KEY)
        except ValueError:
            return None
        return f"{quantity:~.4g}"

    def _selected_units(self) -> tuple[str | None, str | None]:
        if self.data is None:
            return None, None
        x_unit = None
        y_unit = None
        try:
            x_unit = self.data.domain[self.x_variable].attributes.get("unit")
        except KeyError:
            pass
        try:
            y_unit = self.data.domain[self.y_variable].attributes.get("unit")
        except KeyError:
            pass
        return x_unit, y_unit

    def _build_plot(self) -> None:
        self.plot = pg.PlotWidget(self.mainArea)
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        self.mainArea.layout().addWidget(self.plot)

    @Inputs.data
    def set_data(self, data: Table | None) -> None:
        self.data = data
        self.Error.clear()
        self._fit_curve = None
        self.Outputs.fit_results.send(None)
        if self.data is not None:
            self._update_variable_controls()
        self._redraw()

    def _update_variable_controls(self) -> None:
        names = (
            []
            if self.data is None
            else [
                variable.name
                for variable in self.data.domain.attributes
                if isinstance(variable, ContinuousVariable)
            ]
        )
        self.variable_names = names

        x_selected = self._preferred_variable(
            self.x_variable, (DEFAULT_X, FALLBACK_X), names
        )
        y_selected = self._preferred_variable(
            self.y_variable, (DEFAULT_Y, FALLBACK_Y), names
        )
        if y_selected == x_selected and len(names) > 1:
            y_selected = next(name for name in names if name != x_selected)

        self._updating_controls = True
        for combo, selected in (
            (self.x_combo, x_selected),
            (self.y_combo, y_selected),
        ):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(names)
            if selected:
                combo.setCurrentText(selected)
            combo.blockSignals(False)
        # Assigning to an empty combo makes Orange warn; keep the old value.
        if x_selected:
            self.x_variable = x_selected
        if y_selected:
            self.y_variable = y_selected
        self._updating_controls = False

        if self.data is not None and not names:
            self.Error.no_continuous_data()

    @staticmethod
    def _preferred_variable(
        current: str, defaults: tuple[str, ...], names: list[str]
    ) -> str:
        if current in names:
            return current
        for default in defaults:
            if default in names:
                return default
        return names[0] if names else ""

    def _selection_changed(self) -> None:
        if self._updating_controls:
            return
        self._fit_curve = None
        self.Outputs.fit_results.send(None)
        self.Error.fit_failed.clear()
        self._redraw()

    def _xy_data(self) -> tuple[np.ndarray, np.ndarray]:
        if self.data is None:
            raise ValueError("No input data")

        try:
            x_variable = self.data.domain[self.x_variable]
            y_variable = self.data.domain[self.y_variable]
        except KeyError as exc:
            raise ValueError("Select valid x and y variables") from exc

        if not isinstance(x_variable, ContinuousVariable) or not isinstance(
            y_variable, ContinuousVariable
        ):
            # ValueError on purpose: callers treat every bad selection alike.
            raise ValueError("The selected x and y variables must be continuous")  # noqa: TRY004

        _x = np.asarray(self.data.get_column(x_variable), dtype=float)
        _y = np.asarray(self.data.get_column(y_variable), dtype=float)

        # important for lmfit to ensure the values are sorted beforehand
        return (
            _x[np.argsort(_x)],
            _y[
                np.argsort(_x)
            ],  # *1e6 # TODO: clear up this unit business, struggles to fit at natural intensity
        )

    def _redraw(self) -> None:
        self.plot.clear()
        x_unit, y_unit = self._selected_units()
        self.plot.setLabel("bottom", self.x_variable, units=self._unit_symbol(x_unit))
        y_label = self.y_variable
        wavelength = self._measurement_wavelength()
        if wavelength and y_label:
            y_label = f"{y_label} at {wavelength}"
        self.plot.setLabel("left", y_label, units=self._unit_symbol(y_unit))
        if self.data is None or not self.x_variable or not self.y_variable:
            return
        try:
            x, y = self._xy_data()
        except ValueError:
            return
        valid = np.isfinite(x) & np.isfinite(y)
        self.plot.plot(
            x[valid],
            y[valid],
            pen=None,
            symbol="o",
            symbolSize=self.point_size,
            symbolBrush=pg.mkBrush("#1f77b4"),
            symbolPen=pg.mkPen("#1f77b4"),
        )
        if self._fit_curve is not None:
            fit_x = self._fit_curve["x_plt"]
            fit_y = self._fit_curve["y_plt"]
            y_err = self._fit_curve["y_err"]

            # Best-fit line
            self.plot.plot(
                fit_x,
                fit_y,
                pen=pg.mkPen("#d62728", width=self.line_width),
            )

            # Upper/lower bounds of confidence interval
            upper_curve = pg.PlotCurveItem(
                fit_x,
                fit_y + y_err,
                pen=None,
            )

            lower_curve = pg.PlotCurveItem(
                fit_x,
                fit_y - y_err,
                pen=None,
            )

            self.plot.addItem(upper_curve)
            self.plot.addItem(lower_curve)

            # Filled uncertainty region
            fill = pg.FillBetweenItem(
                upper_curve,
                lower_curve,
                brush=pg.mkBrush(214, 39, 40, 80),  # RGBA, alpha ~30%
            )

            self.plot.addItem(fill)
        self.plot.enableAutoRange()

    def clear_fit(self) -> None:
        self._fit_curve = None
        self.Error.fit_failed.clear()
        self.Outputs.fit_results.send(None)
        self._redraw()

    @property
    def binding_model(self) -> BindingModel:
        """The selected binding model."""
        return MODELS[self.fit_model]

    def _select_fit_model_changed(self) -> None:
        self.equation.setText(self.binding_model.equation)

    def fit(self) -> None:
        self.Error.clear()
        if self.data is None:
            self.Outputs.fit_results.send(None)
            return
        try:
            x, y = self._xy_data()
            x_unit, y_unit = self._selected_units()
            fit_result, self._fit_curve = self.binding_model.fit(
                x, y, x_unit=x_unit, y_unit=y_unit
            )
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            self._fit_curve = None
            self.Error.fit_failed(str(exc))
            self.Outputs.fit_results.send(None)
            self._redraw()
            return

        self._redraw()

        self.Outputs.fit_results.send(fit_result)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWBindingPlot).run()
