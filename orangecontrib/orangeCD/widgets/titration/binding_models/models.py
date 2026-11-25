import numpy as np
from lmfit import Model
from lmfit.model import ModelResult
from Orange.data import ContinuousVariable, Domain, StringVariable, Table

# Fitted parameters that share the physical unit of the response (y) or of
# the independent variable (x); everything else (Hill coefficients, r_sq,
# red_chi_sq) is dimensionless.
Y_UNIT_PARAMETERS = {"v_max", "top", "bottom", "p_m"}
X_UNIT_PARAMETERS = {"half_saturation", "k_a", "k_i"}

class BaseModel:

    def __init__(self, x: np.ndarray, y: np.ndarray,
                 model_func: function, eqn_str: str) -> None:
        self.x = x
        self.y = y
        self.model_func = model_func
        self.equation = eqn_str

    def _validate_data(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:

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

        return (np.asarray(x, dtype=float), np.asarray(y, dtype=float))


    def _parameter_unit(self, name: str, x_unit: str | None, y_unit: str | None) -> str:
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
    ):

        parameter_names = list(lmfit_result.params.keys()) + ["r_sq", "red_chi_sq"]
        vals = [i.value for i in lmfit_result.params.values()] + [
            lmfit_result.rsquared,
            lmfit_result.redchi,
        ]
        errs = [i.stderr for i in lmfit_result.params.values()] + [np.nan, np.nan]
        units = [self._parameter_unit(name, x_unit, y_unit) for name in parameter_names]

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
            metas=np.column_stack(
                (
                    np.asarray(parameter_names, dtype=object),
                    np.asarray(units, dtype=object),
                )
            ),
        )
        result.name = "Curve fit"

        return result

    def fitter(
            self,
        x: np.ndarray,
        y: np.ndarray,
        x_unit: str | None = None,
        y_unit: str | None = None,
    ):

        x, y = self._validate_data(x, y)

        top_guess = float(y[np.argmax(x)])
        half_guess = float(np.median(x))

        model = Model(hill_equation)
        params = model.make_params(
            v_max=top_guess,
            half_saturation={"value": half_guess, "min": 0},
            hill_coefficient={"value": 1.0, "min": 0},
        )
        result = model.fit(y, params=params, x=x)

        x_plt = np.linspace(x[0], x[-1])
        y_plt = result.eval(x=x_plt)
        y_err = result.eval_uncertainty(x=x_plt)

        return self.fit_result(result, x_unit, y_unit), {
            "x_plt": x_plt,
            "y_plt": y_plt,
            "y_err": y_err,
        }

    # def obj_func(self, function) -> 



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
    params = model.make_params(
        v_max=top_guess,
        half_saturation={"value": half_guess, "min": 0},
        hill_coefficient={"value": 1.0, "min": 0},
    )
    result = model.fit(y, params=params, x=x)

    x_plt = np.linspace(x[0], x[-1])
    y_plt = result.eval(x=x_plt)
    y_err = result.eval_uncertainty(x=x_plt)

    return fit_result(result, x_unit, y_unit), {
        "x_plt": x_plt,
        "y_plt": y_plt,
        "y_err": y_err,
    }


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
    params = model.make_params(
        bottom=bottom_guess,
        top=top_guess,
        half_saturation={"value": half_guess, "min": 0},
        hill_coefficient={"value": 1.0, "min": 0},
    )
    result = model.fit(y, params=params, x=x)

    x_plt = np.linspace(x[0], x[-1])
    y_plt = result.eval(x=x_plt)
    y_err = result.eval_uncertainty(x=x_plt)

    return fit_result(result, x_unit, y_unit), {
        "x_plt": x_plt,
        "y_plt": y_plt,
        "y_err": y_err,
    }


def fit_bihill(
    x: np.ndarray,
    y: np.ndarray,
    x_unit: str | None = None,
    y_unit: str | None = None,
):

    x, y = validate_data(x, y)

    model = Model(bihill_equation)
    params = model.make_params(
        p_m={"value": y.max(), "min": 0},
        k_a={"value": x[np.argmax(np.gradient(y))], "min": 0},
        k_i={"value": x[np.argmin(np.gradient(y))], "min": 0},
        # possibly could estimate this better in future?
        h_a={
            "value": 1,
            "min": 0,
        },
        h_i={
            "value": 1,
            "min": 0,
        },
    )
    result = model.fit(y, params=params, x=x)

    x_plt = np.linspace(x[0], x[-1])
    y_plt = result.eval(x=x_plt)
    y_err = result.eval_uncertainty(x=x_plt)

    return fit_result(result, x_unit, y_unit), {
        "x_plt": x_plt,
        "y_plt": y_plt,
        "y_err": y_err,
    }


# def validate_data(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:

#     x = np.asarray(x, dtype=float)
#     y = np.asarray(y, dtype=float)
#     valid = np.isfinite(x) & np.isfinite(y)
#     x = x[valid]
#     y = y[valid]

#     if x.size < 4:
#         raise ValueError("At least four finite data points are required for a Hill fit")
#     if np.any(x < 0):
#         raise ValueError("Hill fitting requires non-negative x values")
#     if np.unique(x).size < 4:
#         raise ValueError("At least four distinct x values are required for a Hill fit")

#     positive_x = x[x > 0]
#     if positive_x.size == 0:
#         raise ValueError("At least one x value must be greater than zero")

#     return (x, y)


# def _parameter_unit(name: str, x_unit: str | None, y_unit: str | None) -> str:
#     if name in Y_UNIT_PARAMETERS:
#         return y_unit or ""
#     if name in X_UNIT_PARAMETERS:
#         return x_unit or ""
#     return ""


# def fit_result(
#     lmfit_result: ModelResult,
#     x_unit: str | None = None,
#     y_unit: str | None = None,
# ):

#     parameter_names = list(lmfit_result.params.keys()) + ["r_sq", "red_chi_sq"]
#     vals = [i.value for i in lmfit_result.params.values()] + [
#         lmfit_result.rsquared,
#         lmfit_result.redchi,
#     ]
#     errs = [i.stderr for i in lmfit_result.params.values()] + [np.nan, np.nan]
#     units = [_parameter_unit(name, x_unit, y_unit) for name in parameter_names]

#     domain = Domain(
#         [
#             ContinuousVariable("Estimate"),
#             ContinuousVariable("Standard Error"),
#         ],
#         metas=[StringVariable("Parameter"), StringVariable("Unit")],
#     )
#     result = Table.from_numpy(
#         domain,
#         np.column_stack((vals, errs)),
#         metas=np.column_stack(
#             (
#                 np.asarray(parameter_names, dtype=object),
#                 np.asarray(units, dtype=object),
#             )
#         ),
#     )
#     result.name = "Curve fit"

#     return result

