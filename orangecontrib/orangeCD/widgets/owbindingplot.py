#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Scatter plot and manual Hill-equation fitting for Binding Data tables."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from scipy.optimize import curve_fit

from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget


DEFAULT_X = "Titration point"
DEFAULT_Y = "Delta A"


def hill_equation(
    x: np.ndarray,
    bottom: float,
    top: float,
    half_saturation: float,
    hill_coefficient: float,
) -> np.ndarray:
    """Four-parameter Hill equation.

    y = bottom + (top - bottom) * x**n / (K_half**n + x**n)
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        x_power = np.power(x, hill_coefficient)
        half_power = np.power(half_saturation, hill_coefficient)
        return bottom + (top - bottom) * x_power / (half_power + x_power)


def fit_hill_equation(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Fit a four-parameter Hill model and return parameters and diagnostics."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if x.size < 4:
        raise ValueError("At least four finite data points are required for a Hill fit")
    if np.any(x < 0):
        raise ValueError("Hill fitting requires non-negative x values")
    if np.unique(x).size < 4:
        raise ValueError("At least four distinct x values are required for a Hill fit")

    positive_x = x[x > 0]
    if positive_x.size == 0:
        raise ValueError("At least one x value must be greater than zero")

    bottom_guess = float(y[np.argmin(x)])
    top_guess = float(y[np.argmax(x)])
    half_guess = float(np.median(positive_x))
    span = max(float(np.ptp(y)), abs(top_guess), abs(bottom_guess), 1.0)

    parameters, covariance = curve_fit(
        hill_equation,
        x,
        y,
        p0=(bottom_guess, top_guess, half_guess, 1.0),
        bounds=(
            (-np.inf, -np.inf, np.finfo(float).eps, 0.01),
            (np.inf, np.inf, np.inf, 20.0),
        ),
        maxfev=50000,
    )

    fitted = hill_equation(x, *parameters)
    residuals = y - fitted
    residual_sum_squares = float(np.sum(residuals**2))
    total_sum_squares = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = (
        1.0 - residual_sum_squares / total_sum_squares
        if total_sum_squares > 0
        else np.nan
    )
    rmse = float(np.sqrt(np.mean(residuals**2)))
    return parameters, covariance, r_squared, rmse


class OWBindingPlot(OWWidget):
    name = "Binding Plot"
    description = "Plot Binding Data and manually fit a four-parameter Hill equation."
    icon = "icons/Titration.svg"
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

    class Error(OWWidget.Error):
        no_continuous_data = Msg("Input contains no continuous variables.")
        invalid_selection = Msg("Select valid x and y variables.")
        fit_failed = Msg("{}")

    def __init__(self) -> None:
        super().__init__()
        self.data: Table | None = None
        self.variable_names: list[str] = []
        self._updating_controls = False
        self._fit_curve: tuple[np.ndarray, np.ndarray] | None = None
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
            valueType=str,
            orientation="horizontal",
            callback=self._selection_changed,
        )
        self.y_combo = gui.comboBox(
            variables_box,
            self,
            "y_variable",
            label="Y variable",
            items=[],
            sendSelectedValue=True,
            valueType=str,
            orientation="horizontal",
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
            orientation="horizontal",
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
            orientation="horizontal",
            callback=self._redraw,
        )

        fit_box = gui.widgetBox(self.controlArea, "Hill fit")
        equation = gui.widgetLabel(
            fit_box,
            "y = bottom + (top - bottom) x^n / (K_half^n + x^n)",
        )
        equation.setWordWrap(True)
        gui.button(
            fit_box,
            self,
            "Fit Hill equation",
            callback=self.fit,
        )
        gui.button(
            fit_box,
            self,
            "Clear fit",
            callback=self.clear_fit,
        )
        gui.rubber(self.controlArea)

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
        self._update_variable_controls()
        self._redraw()

    def _update_variable_controls(self) -> None:
        names = [] if self.data is None else [
            variable.name
            for variable in self.data.domain.attributes
            if isinstance(variable, ContinuousVariable)
        ]
        self.variable_names = names

        x_selected = self._preferred_variable(self.x_variable, DEFAULT_X, names)
        y_selected = self._preferred_variable(self.y_variable, DEFAULT_Y, names)
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
        self.x_variable = x_selected
        self.y_variable = y_selected
        self._updating_controls = False

        if self.data is not None and not names:
            self.Error.no_continuous_data()

    @staticmethod
    def _preferred_variable(current: str, default: str, names: list[str]) -> str:
        if current in names:
            return current
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
            raise ValueError("The selected x and y variables must be continuous")
        return (
            np.asarray(self.data.get_column(x_variable), dtype=float),
            np.asarray(self.data.get_column(y_variable), dtype=float),
        )

    def _redraw(self) -> None:
        self.plot.clear()
        self.plot.setLabel("bottom", self.x_variable)
        self.plot.setLabel("left", self.y_variable)
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
            fit_x, fit_y = self._fit_curve
            self.plot.plot(
                fit_x,
                fit_y,
                pen=pg.mkPen("#d62728", width=self.line_width),
            )
        self.plot.enableAutoRange()

    def clear_fit(self) -> None:
        self._fit_curve = None
        self.Error.fit_failed.clear()
        self.Outputs.fit_results.send(None)
        self._redraw()

    def fit(self) -> None:
        self.Error.clear()
        if self.data is None:
            self.Outputs.fit_results.send(None)
            return
        try:
            x, y = self._xy_data()
            parameters, covariance, r_squared, rmse = fit_hill_equation(x, y)
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            self._fit_curve = None
            self.Error.fit_failed(str(exc))
            self.Outputs.fit_results.send(None)
            self._redraw()
            return

        valid_x = x[np.isfinite(x) & np.isfinite(y)]
        fit_x = np.linspace(float(np.min(valid_x)), float(np.max(valid_x)), 500)
        self._fit_curve = (fit_x, hill_equation(fit_x, *parameters))
        self._redraw()

        standard_errors = np.sqrt(np.diag(covariance))
        parameter_names = [
            "bottom",
            "top",
            "K_half",
            "Hill coefficient",
            "R squared",
            "RMSE",
        ]
        estimates = np.concatenate((parameters, [r_squared, rmse]))
        errors = np.concatenate((standard_errors, [np.nan, np.nan]))

        domain = Domain(
            [
                ContinuousVariable("Estimate"),
                ContinuousVariable("Standard Error"),
            ],
            metas=[StringVariable("Parameter")],
        )
        result = Table.from_numpy(
            domain,
            np.column_stack((estimates, errors)),
            metas=np.asarray(parameter_names, dtype=object).reshape(-1, 1),
        )
        result.name = f"Hill fit: {self.y_variable} vs {self.x_variable}"
        result.attributes.update(
            {
                "model": (
                    "bottom + (top - bottom) * x^n / (K_half^n + x^n)"
                ),
                "x_variable": self.x_variable,
                "y_variable": self.y_variable,
                "point_count": int(valid_x.size),
            }
        )
        self.Outputs.fit_results.send(result)


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWBindingPlot).run()
