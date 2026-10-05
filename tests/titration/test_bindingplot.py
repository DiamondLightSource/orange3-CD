import numpy as np
import pytest
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets.tests import base as orange_tests

from orangecontrib.orangeCD.widgets.titration.owbindingplot import (
    OWBindingPlot,
    bihill_equation,
    fit_bihill,
    fit_hill,
    fit_hill1,
    fit_result,
    hill1_equation,
    hill_equation,
    validate_data,
)
from orangecontrib.orangeCD.widgets.titration.utils import (
    MEASUREMENT_WAVELENGTH_KEY,
)

X = np.array([0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0])


def parameters(table):
    return dict(zip(table.metas[:, 0], table.X[:, 0])), dict(
        zip(table.metas[:, 0], table.metas[:, 1])
    )


class TestEquations:
    def test_hill_half_saturation(self):
        assert hill_equation(np.array([5.0]), 10.0, 5.0, 2.0)[0] == pytest.approx(5.0)

    def test_hill_limits(self):
        y = hill_equation(np.array([1e-6, 1e6]), 10.0, 5.0, 1.0)
        assert y[0] == pytest.approx(0, abs=1e-5)
        assert y[1] == pytest.approx(10.0, rel=1e-4)

    def test_hill1_bounds(self):
        y = hill1_equation(np.array([1e-9, 1e9]), 2.0, 8.0, 5.0, 1.0)
        assert y[0] == pytest.approx(2.0)
        assert y[1] == pytest.approx(8.0)
        assert hill1_equation(np.array([5.0]), 2.0, 8.0, 5.0, 3.0)[0] == pytest.approx(
            5.0
        )

    def test_bihill_peaks_between_constants(self):
        x = np.logspace(-2, 3, 200)
        y = bihill_equation(x, 10.0, 1.0, 2.0, 100.0, 2.0)
        assert 1.0 < x[np.argmax(y)] < 100.0
        assert y.max() <= 10.0

    def test_equations_accept_zero_without_raising(self):
        with np.errstate(all="raise"):
            hill_equation([0.0], 1, 1, 1)
            bihill_equation([0.0], 1, 1, 1, 1, 1)


class TestValidateData:
    def test_drops_non_finite(self):
        x, y = validate_data([1, 2, 3, 4, 5, np.nan], [1, 2, 3, 4, np.inf, 6])
        assert list(x) == [1, 2, 3, 4]
        assert list(y) == [1, 2, 3, 4]

    @pytest.mark.parametrize(
        "x, message",
        [
            ([1, 2, 3], "four"),
            ([1, 2, 2, 2, 2], "distinct"),
            ([-1, 1, 2, 3], "non-negative"),
        ],
    )
    def test_invalid(self, x, message):
        with pytest.raises(ValueError, match=message):
            validate_data(x, np.ones(len(x)))


class TestFitting:
    def test_hill_recovers_parameters(self):
        y = hill_equation(X, 10.0, 2.0, 1.5)
        result, curve = fit_hill(X, y, x_unit="molar", y_unit="millidegree")
        values, units = parameters(result)
        assert values["v_max"] == pytest.approx(10.0, rel=1e-3)
        assert values["half_saturation"] == pytest.approx(2.0, rel=1e-3)
        assert values["hill_coefficient"] == pytest.approx(1.5, rel=1e-3)
        assert values["r_sq"] == pytest.approx(1.0)
        assert units["v_max"] == "millidegree"
        assert units["half_saturation"] == "molar"
        assert units["hill_coefficient"] == ""
        assert set(curve) == {"x_plt", "y_plt", "y_err"}
        assert len(curve["x_plt"]) == len(curve["y_plt"]) == len(curve["y_err"])

    def test_hill1_recovers_parameters(self):
        y = hill1_equation(X, 1.0, 9.0, 2.0, 1.0)
        values, units = parameters(fit_hill1(X, y, "molar", "mdeg")[0])
        assert values["bottom"] == pytest.approx(1.0, rel=1e-2)
        assert values["top"] == pytest.approx(9.0, rel=1e-2)
        assert units["bottom"] == units["top"] == "mdeg"

    def test_bihill_fits(self):
        x = np.logspace(-1, 2, 25)
        y = bihill_equation(x, 10.0, 1.0, 2.0, 30.0, 2.0)
        values, units = parameters(fit_bihill(x, y, "molar", "mdeg")[0])
        assert values["r_sq"] > 0.99
        assert units["p_m"] == "mdeg"
        assert units["k_a"] == units["k_i"] == "molar"

    def test_unsorted_inputs_still_fit_after_sorting(self):
        order = np.argsort(-X)
        y = hill_equation(X, 4.0, 3.0, 1.0)
        result, _ = fit_hill(X[order][::-1], y[order][::-1])
        assert parameters(result)[0]["v_max"] == pytest.approx(4.0, rel=1e-3)

    def test_invalid_data_raises(self):
        with pytest.raises(ValueError):
            fit_hill([1, 2], [1, 2])

    def test_fit_result_table_layout(self):
        result, _ = fit_hill(X, hill_equation(X, 10.0, 2.0, 1.0))
        assert result.name == "Curve fit"
        assert [v.name for v in result.domain.attributes] == [
            "Estimate",
            "Standard Error",
        ]
        assert [v.name for v in result.domain.metas] == ["Parameter", "Unit"]
        assert list(result.metas[:, 0])[-2:] == ["r_sq", "red_chi_sq"]

    def test_fit_result_function_exposed(self):
        assert callable(fit_result)


def binding_table(y=None, wavelength=None, n=7, with_units=True):
    titration_point = np.linspace(0, 3, n)
    y = hill_equation(titration_point, 10.0, 1.0, 1.0) if y is None else y
    domain = Domain(
        [
            ContinuousVariable("Titration point"),
            ContinuousVariable("Delta A"),
            ContinuousVariable("CD"),
        ],
        metas=[StringVariable("Sample")],
    )
    if with_units:
        domain["Delta A"].attributes["unit"] = "millidegree"
        domain["Titration point"].attributes["unit"] = "molar"
    table = Table.from_numpy(
        domain,
        np.column_stack((titration_point, y, np.arange(n, dtype=float))),
        metas=np.array([[f"s{i}"] for i in range(n)], dtype=object),
    )
    if wavelength is not None:
        table.attributes[MEASUREMENT_WAVELENGTH_KEY] = wavelength
    return table


class TestWidget(orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWBindingPlot)

    def send(self, table):
        self.send_signal(self.widget.Inputs.data, table, widget=self.widget)

    def output(self):
        return self.get_output(self.widget.Outputs.fit_results, widget=self.widget)

    def test_default_variables(self):
        self.send(binding_table())
        assert self.widget.x_variable == "Titration point"
        assert self.widget.y_variable == "Delta A"
        assert self.widget.variable_names == ["Titration point", "Delta A", "CD"]

    def test_y_differs_from_x_when_defaults_missing(self):
        domain = Domain([ContinuousVariable("p"), ContinuousVariable("q")])
        self.send(Table.from_numpy(domain, np.ones((5, 2))))
        assert self.widget.x_variable == "p" and self.widget.y_variable == "q"

    def test_no_continuous_data(self):
        table = Table.from_numpy(
            Domain([], metas=[StringVariable("m")]),
            np.zeros((3, 0)),
            metas=np.array([["a"], ["b"], ["c"]], dtype=object),
        )
        self.send(table)
        assert self.widget.Error.no_continuous_data.is_shown()

    def test_points_are_plotted_sorted(self):
        table = binding_table()
        self.send(table)
        (item,) = self.widget.plot.listDataItems()
        x, _ = item.getData()
        assert list(x) == sorted(x)
        assert len(x) == 7

    def test_axis_labels_with_units(self):
        self.send(binding_table(wavelength="260.0 nanometer"))
        assert self.widget.plot.getAxis("bottom").labelUnits == "M"
        left = self.widget.plot.getAxis("left")
        assert left.labelUnits == "mdeg"
        assert left.labelText == "Delta A at 260 nm"

    def test_measurement_wavelength_helper(self):
        self.send(binding_table())
        assert self.widget._measurement_wavelength() is None
        self.send(binding_table(wavelength="garbage"))
        assert self.widget._measurement_wavelength() is None
        self.send(binding_table(wavelength="402.5 nanometer"))
        assert self.widget._measurement_wavelength() == "402.5 nm"
        self.send(None)
        assert self.widget._measurement_wavelength() is None

    def test_fit_sends_results_and_draws_curve(self):
        self.send(binding_table())
        self.widget.fit()
        out = self.output()
        assert out is not None and out.name == "Curve fit"
        values, units = parameters(out)
        assert values["v_max"] == pytest.approx(10.0, rel=1e-2)
        assert units["v_max"] == "millidegree"
        assert units["half_saturation"] == "molar"
        assert self.widget._fit_curve is not None
        assert len(self.widget.plot.listDataItems()) >= 2

    def test_fit_failure_is_reported(self):
        self.send(binding_table(n=3))
        self.widget.fit()
        assert self.widget.Error.fit_failed.is_shown()
        assert self.output() is None
        assert self.widget._fit_curve is None

    def test_fit_without_data(self):
        self.widget.fit()
        assert self.output() is None
        assert not self.widget.Error.active

    def test_clear_fit(self):
        self.send(binding_table())
        self.widget.fit()
        self.widget.clear_fit()
        assert self.widget._fit_curve is None
        assert self.output() is None

    def test_changing_variables_discards_fit(self):
        self.send(binding_table())
        self.widget.fit()
        self.widget.y_variable = "CD"
        self.widget._selection_changed()
        assert self.widget._fit_curve is None
        assert self.output() is None

    def test_new_data_discards_fit(self):
        self.send(binding_table())
        self.widget.fit()
        self.send(binding_table())
        assert self.widget._fit_curve is None
        assert self.output() is None

    def test_model_switch(self):
        self.send(binding_table())
        self.widget.fit_model = 1
        self.widget._select_fit_model_changed()
        assert self.widget.fit_model_fitter is fit_hill1
        assert "bottom" in self.widget.equation.text()
        self.widget.fit()
        values, _ = parameters(self.output())
        assert "bottom" in values and "top" in values

    def test_non_continuous_selection_is_invalid(self):
        self.send(binding_table())
        with self.assertWarns(UserWarning):  # Orange warns about the stale value
            self.widget.x_variable = "Sample"
        with pytest.raises(ValueError):
            self.widget._xy_data()
