#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Scatter plot and manual Hill-equation fitting for Binding Data tables."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from scipy.optimize import curve_fit
from lmfit import Model
from lmfit.model import ModelResult

from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from . import Q_

DEFAULT_X = "Titration point"
DEFAULT_Y = "Delta A"

# Fitted parameters that share the physical unit of the response (y) or of
# the independent variable (x); everything else (Hill coefficients, r_sq,
# red_chi_sq) is dimensionless.
Y_UNIT_PARAMETERS = {"v_max", "top", "bottom", "p_m"}
X_UNIT_PARAMETERS = {"half_saturation", "k_a", "k_i"}

"""
functions for fitting curves. These are pure function models, of form
func(x: np.adarray, *args) -> nd.array, where func() performs some operation
on x with the values in args
"""
def hill_equation(
    x: np.ndarray,
    v_max: float,
    half_saturation: float,
    hill_coefficient: float,
) -> np.ndarray:
    """Three-parameter Hill equation.

    y = v_max * (x**n / (K_half**n + x**n))
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        x_power = np.power(x, hill_coefficient)
        half_power = np.power(half_saturation, hill_coefficient)
        return v_max * x_power / (half_power + x_power)

def hill1_equation(
    x: np.ndarray,
    bottom: float,
    top: float,
    half_saturation: float,
    hill_coefficient: float,
) -> np.ndarray:
    """Four-parameter Hill1 equation.

    y = bottom + (top - bottom) * (x**n / (K_half**n + x**n))
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        x_power = np.power(x, hill_coefficient)
        half_power = np.power(half_saturation, hill_coefficient)
        return bottom + (top - bottom) * x_power / (half_power + x_power)

def bihill_equation(
    x: np.ndarray,
    p_m: float,
    k_a: float,
    h_a: float,
    k_i: float,
    h_i: float
) -> np.ndarray:
    """Five-parameter BiHill equation

    y = p_m / ((1 + ((k_a/x)^h_a))*(1 + ((x/k_i)^h_i)))
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        term1 = 1 + np.power((k_a / x), h_a)
        term2 = 1 + np.power((x / k_i), h_i)
        den = term1 * term2
        return p_m / den

"""
Fitting wrapper functions. These take the data to fit, and perform the fitting
using lmfit, constructing the model and parameters using appropriate estimates.
Before fitting, they check the input data is valid using validate_data().
After fitting, a Table is constructed using the fit_result() function.
"""
def fit_hill(
        x: np.ndarray,
        y: np.ndarray,
        x_unit: str | None = None,
        y_unit: str | None = None,
        ):

    x, y = validate_data(x, y)

    top_guess = float(y[np.argmax(x)])
    half_guess = float(np.median(x))

    model = Model(hill_equation)
    params = model.make_params(v_max = top_guess,
                               half_saturation = dict(value = half_guess, min = 0),
                               hill_coefficient = dict(value = 1.0, min = 0)
                               )
    result = model.fit(y, params=params, x=x)

    x_plt = np.linspace(x[0], x[-1])
    y_plt = result.eval(x=x_plt)
    y_err = result.eval_uncertainty(x=x_plt)

    return fit_result(result, x_unit, y_unit), {'x_plt': x_plt, 'y_plt': y_plt, 'y_err': y_err}

def fit_hill1(
        x: np.ndarray,
        y: np.ndarray,
        x_unit: str | None = None,
        y_unit: str | None = None,
        ):
    x, y = validate_data(x, y)

    bottom_guess = float(y[np.argmin(x)])
    top_guess = float(y[np.argmax(x)])
    half_guess = float(np.median(x))

    model = Model(hill1_equation)
    params = model.make_params(bottom = bottom_guess,
                               top = top_guess,
                               half_saturation = dict(value = half_guess, min = 0),
                               hill_coefficient = dict(value = 1.0, min = 0)
                               )
    result = model.fit(y, params=params, x=x)

    x_plt = np.linspace(x[0], x[-1])
    y_plt = result.eval(x=x_plt)
    y_err = result.eval_uncertainty(x=x_plt)

    return fit_result(result, x_unit, y_unit), {'x_plt': x_plt, 'y_plt': y_plt, 'y_err': y_err}

def fit_bihill(
        x: np.ndarray,
        y: np.ndarray,
        x_unit: str | None = None,
        y_unit: str | None = None,
        ):

    x, y = validate_data(x, y)

    model = Model(bihill_equation)
    params = model.make_params(p_m = dict(value = y.max(), min = 0),
                               k_a = dict(value = x[np.argmax(np.gradient(y))], min = 0),
                               k_i = dict(value = x[np.argmin(np.gradient(y))], min = 0),
                               # possibly could estimate this better in future?
                               h_a = dict(value = 1, min = 0,),
                               h_i = dict(value = 1, min = 0,),
                               )
    result = model.fit(y, params=params, x=x)

    x_plt = np.linspace(x[0], x[-1])
    y_plt = result.eval(x=x_plt)
    y_err = result.eval_uncertainty(x=x_plt)

    return fit_result(result, x_unit, y_unit), {'x_plt': x_plt, 'y_plt': y_plt, 'y_err': y_err}

def validate_data(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:

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

    return (x, y)

def _parameter_unit(name: str, x_unit: str | None, y_unit: str | None) -> str:
    if name in Y_UNIT_PARAMETERS:
        return y_unit or ""
    if name in X_UNIT_PARAMETERS:
        return x_unit or ""
    return ""


def fit_result(
    lmfit_result: ModelResult,
    x_unit: str | None = None,
    y_unit: str | None = None,
):

    parameter_names = list(lmfit_result.params.keys()) + ['r_sq', 'red_chi_sq']
    vals = [i.value for i in lmfit_result.params.values()] + [lmfit_result.rsquared, lmfit_result.redchi]
    errs = [i.stderr for i in lmfit_result.params.values()] + [np.nan, np.nan]
    units = [_parameter_unit(name, x_unit, y_unit) for name in parameter_names]

    domain = Domain(
        [
            ContinuousVariable("Estimate"),
            ContinuousVariable("Standard Error"),
        ],
        metas=[StringVariable("Parameter"), StringVariable("Unit")],
    )
    result = Table.from_numpy(
        domain,
        np.column_stack((vals, errs)),
        metas=np.column_stack((
            np.asarray(parameter_names, dtype=object),
            np.asarray(units, dtype=object),
        )),
    )
    result.name = "Curve fit"

    return result


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
    fitter_function = Setting(0)
    fit_model = Setting("Hill model (3 variable)")


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
        self.fit_model_dict = {
            0: {
                "text": "Hill model (3 variable)",
                "function": hill_equation, 
                "fitter": fit_hill, 
                "equation": "y = v_max * (x**n / (K_half**n + x**n))"
                },
            1: {
                "text": "Hill model (4 variable)",
                "function": hill1_equation,
                "fitter": fit_hill1,
                "equation": "y = bottom + (top - bottom) * (x**n / (K_half**n + x**n))"
                },
            2: {    
                "text": "BiHill model",
                "function": bihill_equation,
                "fitter": fit_bihill,
                "equation": "y = p_m / ((1 + ((k_a/x)^h_a))*(1 + ((x/k_i)^h_i)))"
                }
                }
        self.fit_model_list = tuple(i["text"] for i in self.fit_model_dict.values())
        # need to initialise these here otherwise get mixed up in class, possibly
        self.fit_model_func = self.fit_model_dict.get(self.fitter_function).get("function")
        self.fit_model_fitter = self.fit_model_dict.get(self.fitter_function).get("fitter")
        self.equation_str = self.fit_model_dict.get(self.fitter_function).get("equation")

        self.equation = gui.label(
            fit_box,
            self,
            f"{self.equation_str}"
        )
        self.equation.setWordWrap(True)

        gui.comboBox(
            fit_box,
            self,
            "fit_model",
            label="Select fit model",
            items=self.fit_model_list,
            sendSelectedValue=False,
            valueType=int,
            orientation="horizontal",
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
            y_variable, ContinuousVariable):
            raise ValueError("The selected x and y variables must be continuous")

        _x = np.asarray(self.data.get_column(x_variable), dtype=float)
        _y = np.asarray(self.data.get_column(y_variable), dtype=float)

        # important for lmfit to ensure the values are sorted beforehand
        return (
            _x[np.argsort(_x)],
            _y[np.argsort(_x)]#*1e6 # TODO: clear up this unit business, struggles to fit at natural intensity
        )

    def _redraw(self) -> None:
        self.plot.clear()
        x_unit, y_unit = self._selected_units()
        self.plot.setLabel("bottom", self.x_variable, units=self._unit_symbol(x_unit))
        self.plot.setLabel("left", self.y_variable, units=self._unit_symbol(y_unit))
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
            fit_curve = self.plot.plot(
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

    def _select_fit_model_changed(self) -> None:
        self.fit_model_func = self.fit_model_dict.get(self.fit_model).get("function")
        self.fit_model_fitter = self.fit_model_dict.get(self.fit_model).get("fitter")
        self.equation.setText(self.fit_model_dict.get(self.fit_model).get("equation"))

    def fit(self) -> None:
        self.Error.clear()
        if self.data is None:
            self.Outputs.fit_results.send(None)
            return
        try:
            x, y = self._xy_data()
            x_unit, y_unit = self._selected_units()
            fit_result, self._fit_curve = self.fit_model_fitter(
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
