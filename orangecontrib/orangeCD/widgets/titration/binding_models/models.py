"""Binding models fitted by the Binding Plot widget.

Each model pairs a pure equation, ``func(x, *parameters) -> y``, with the
initial parameter estimates for fitting it. The shared work (validating the
data, fitting with lmfit and building the result table) is done by
:class:`BindingModel`; a new model only needs to subclass it and set
``name``, ``equation``, ``function`` and :meth:`~BindingModel.initial_params`.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import ClassVar

import numpy as np
from lmfit import Model, Parameters
from lmfit.model import ModelResult
from Orange.data import ContinuousVariable, Domain, StringVariable, Table

# Fitted parameters that share the physical unit of the response (y) or of
# the independent variable (x); everything else (Hill coefficients, r_sq,
# red_chi_sq) is dimensionless.
Y_UNIT_PARAMETERS = {"v_max", "top", "bottom", "p_m"}
X_UNIT_PARAMETERS = {"half_saturation", "k_a", "k_i"}

# Number of points on the plotted fit curve.
CURVE_POINTS = 50


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
    x: np.ndarray, p_m: float, k_a: float, h_a: float, k_i: float, h_i: float
) -> np.ndarray:
    """Five-parameter BiHill equation.

    y = p_m / ((1 + ((k_a/x)^h_a))*(1 + ((x/k_i)^h_i)))
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        term1 = 1 + np.power((k_a / x), h_a)
        term2 = 1 + np.power((x / k_i), h_i)
        return p_m / (term1 * term2)


class BindingModel:
    """A binding equation and how to fit it to titration data."""

    #: Label shown in the widget.
    name: ClassVar[str]
    #: Human-readable form of the equation.
    equation: ClassVar[str]
    #: The pure equation, ``function(x, *parameters) -> y``.
    function: ClassVar[Callable[..., np.ndarray]]
    #: Fewest finite (and distinct) x values a fit needs.
    min_points: ClassVar[int] = 4

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"

    def initial_params(self, model: Model, x: np.ndarray, y: np.ndarray) -> Parameters:
        """Starting values (and bounds) of the parameters for this data."""
        raise NotImplementedError

    def validate_data(
        self, x: np.ndarray, y: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Drop non-finite pairs and check there is enough data to fit."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        x, y = x[valid], y[valid]

        if x.size < self.min_points:
            raise ValueError(
                f"The {self.name} fit needs at least {self.min_points} finite "
                "data points"
            )
        if np.any(x < 0):
            raise ValueError(f"The {self.name} fit needs non-negative x values")
        if np.unique(x).size < self.min_points:
            raise ValueError(
                f"The {self.name} fit needs at least {self.min_points} distinct "
                "x values"
            )
        if not np.any(x > 0):
            raise ValueError("At least one x value must be greater than zero")
        return x, y

    def fit(
        self,
        x: np.ndarray,
        y: np.ndarray,
        x_unit: str | None = None,
        y_unit: str | None = None,
    ) -> tuple[Table, dict[str, np.ndarray]]:
        """Fit the model to ``(x, y)``.

        Returns the table of fitted parameters and the fitted curve, as
        ``{"x_plt", "y_plt", "y_err"}`` arrays spanning the range of x.
        Raises ``ValueError`` for unusable data.
        """
        x, y = self.validate_data(x, y)
        model = Model(type(self).function)
        result = model.fit(y, params=self.initial_params(model, x, y), x=x)

        x_plt = np.linspace(x.min(), x.max(), CURVE_POINTS)
        curve = {
            "x_plt": x_plt,
            "y_plt": result.eval(x=x_plt),
            "y_err": result.eval_uncertainty(x=x_plt),
        }
        return self.fit_result(result, x_unit, y_unit), curve

    @staticmethod
    def parameter_unit(name: str, x_unit: str | None, y_unit: str | None) -> str:
        if name in Y_UNIT_PARAMETERS:
            return y_unit or ""
        if name in X_UNIT_PARAMETERS:
            return x_unit or ""
        return ""

    def fit_result(
        self,
        lmfit_result: ModelResult,
        x_unit: str | None = None,
        y_unit: str | None = None,
    ) -> Table:
        """Table of the fitted parameters with their errors and units."""
        params = lmfit_result.params.values()
        names = [*lmfit_result.params.keys(), "r_sq", "red_chi_sq"]
        values = [p.value for p in params] + [
            lmfit_result.rsquared,
            lmfit_result.redchi,
        ]
        # stderr is None when lmfit cannot estimate the uncertainties.
        errors = [np.nan if p.stderr is None else p.stderr for p in params]
        errors += [np.nan, np.nan]
        units = [self.parameter_unit(name, x_unit, y_unit) for name in names]

        domain = Domain(
            [ContinuousVariable("Estimate"), ContinuousVariable("Standard Error")],
            metas=[StringVariable("Parameter"), StringVariable("Unit")],
        )
        result = Table.from_numpy(
            domain,
            np.column_stack((values, errors)).astype(float),
            metas=np.column_stack(
                (np.asarray(names, dtype=object), np.asarray(units, dtype=object))
            ),
        )
        result.name = "Curve fit"
        return result


class HillModel(BindingModel):
    name = "Hill model (3 variable)"
    equation = "y = v_max * (x**n / (K_half**n + x**n))"
    function = staticmethod(hill_equation)

    def initial_params(self, model, x, y):
        return model.make_params(
            v_max=float(y[np.argmax(x)]),
            half_saturation={"value": float(np.median(x)), "min": 0},
            hill_coefficient={"value": 1.0, "min": 0},
        )


class Hill1Model(BindingModel):
    name = "Hill model (4 variable)"
    equation = "y = bottom + (top - bottom) * (x**n / (K_half**n + x**n))"
    function = staticmethod(hill1_equation)

    def initial_params(self, model, x, y):
        return model.make_params(
            bottom=float(y[np.argmin(x)]),
            top=float(y[np.argmax(x)]),
            half_saturation={"value": float(np.median(x)), "min": 0},
            hill_coefficient={"value": 1.0, "min": 0},
        )


class BiHillModel(BindingModel):
    name = "BiHill model"
    equation = "y = p_m / ((1 + ((k_a/x)^h_a))*(1 + ((x/k_i)^h_i)))"
    function = staticmethod(bihill_equation)

    def initial_params(self, model, x, y):
        # Steepest rise and fall, per point: the points must be in x order.
        order = np.argsort(x)
        x, y = x[order], y[order]
        gradient = np.gradient(y)
        return model.make_params(
            p_m={"value": float(y.max()), "min": 0},
            k_a={"value": float(x[np.argmax(gradient)]), "min": 0},
            k_i={"value": float(x[np.argmin(gradient)]), "min": 0},
            # possibly could estimate these better in future?
            h_a={"value": 1.0, "min": 0},
            h_i={"value": 1.0, "min": 0},
        )


#: The models offered by the Binding Plot widget, in menu order. The widget
#: stores the index into this tuple, so append new models at the end.
MODELS: tuple[BindingModel, ...] = (HillModel(), Hill1Model(), BiHillModel())
