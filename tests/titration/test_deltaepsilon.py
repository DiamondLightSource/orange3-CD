import numpy as np
import pytest
from Orange.data import Table
from Orange.widgets.tests.base import WidgetTest
from pint import DimensionalityError

from helpers import spectra_table, titration_table
from orangecontrib.orangeCD.widgets.titration import Q_
from orangecontrib.orangeCD.widgets.titration.owdeltaepsilon import (
    DELTA_EPSILON_UNIT,
    OWDeltaEpsilon,
    calculate_delta_epsilon,
)
from orangecontrib.orangeCD.widgets.titration.utils import (
    CONCENTRATION_KEY,
    PATHLENGTH_KEY,
    spectrum_names,
    spectrum_units,
    table_quantity,
    unit_string,
)

NAMES = ["Background | sol_A", "a | raw_data", "b | raw_data"]
CD = np.array([[-5.0, -5.0, -5.0, -5.0], [-10.0, -20.0, -30.0, -40.0], [-1.0, 2.0, 3.0, 4.0]])
CONC = 15.0  # uM


def expected(cd_mdeg, conc_um=CONC, pathlength=1.0, mrw=113.0, mw_a=1.0):
    return cd_mdeg / 32980 * mrw / mw_a / (conc_um * 1e-6 * pathlength)


class TestCalculation:
    def args(self, **overrides):
        values = dict(
            cd=Q_(np.array([10.0, -20.0]), "millidegree"),
            concentration=Q_(15, "micromolar"),
            pathlength=Q_(1, "centimeter"),
            mean_residue_molecular_weight=Q_(113, "gram / mole"),
            solution_a_molecular_weight=Q_(1, "gram / mole"),
        )
        values.update(overrides)
        return values

    def test_value_and_unit(self):
        result = calculate_delta_epsilon(**self.args())
        np.testing.assert_allclose(result.magnitude, expected(np.array([10.0, -20.0])))
        assert result.units == Q_(1, DELTA_EPSILON_UNIT).units

    def test_degree_equals_millidegree(self):
        in_degree = calculate_delta_epsilon(**self.args(cd=Q_(np.array([0.01, -0.02]), "degree")))
        np.testing.assert_allclose(in_degree.magnitude, expected(np.array([10.0, -20.0])))

    def test_other_units_are_converted(self):
        result = calculate_delta_epsilon(
            **self.args(
                concentration=Q_(0.015, "millimolar"),
                pathlength=Q_(10, "millimeter"),
                mean_residue_molecular_weight=Q_(0.113, "kilogram / mole"),
            )
        )
        np.testing.assert_allclose(result.magnitude, expected(np.array([10.0, -20.0])))

    @pytest.mark.parametrize(
        "name, value",
        [
            ("concentration", Q_(0, "micromolar")),
            ("pathlength", Q_(-1, "centimeter")),
            ("mean_residue_molecular_weight", Q_(float("nan"), "gram / mole")),
            ("solution_a_molecular_weight", Q_(float("inf") * -1, "gram / mole")),
        ],
    )
    def test_non_positive_parameters(self, name, value):
        with pytest.raises(ValueError, match="greater than zero"):
            calculate_delta_epsilon(**self.args(**{name: value}))

    def test_incompatible_cd_unit(self):
        with pytest.raises(DimensionalityError):
            calculate_delta_epsilon(**self.args(cd=Q_(np.array([1.0]), "nanometer")))


class TestWidget(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWDeltaEpsilon)

    def run_widget(self, table=None, titration=None):
        table = spectra_table(NAMES, CD) if table is None else table
        titration = titration_table([0.5, 1.0]) if titration is None else titration
        self.send_signal(self.widget.Inputs.titration, titration, widget=self.widget)
        self.send_signal(self.widget.Inputs.data, table, widget=self.widget)
        return self.get_output(self.widget.Outputs.data, widget=self.widget)

    def test_output(self):
        out = self.run_widget()
        assert not self.widget.Error.active
        assert spectrum_names(out) == [
            *NAMES,
            "Background | sol_A_delta_epsilon",
            "a | raw_data_delta_epsilon",
            "b | raw_data_delta_epsilon",
        ]
        assert spectrum_units(out) == ["millidegree"] * 3 + [unit_string(DELTA_EPSILON_UNIT)] * 3
        np.testing.assert_allclose(out.X[:3], CD)
        np.testing.assert_allclose(out.X[3:], expected(CD))

    def test_output_attributes(self):
        out = self.run_widget()
        assert out.attributes["wavelength_unit"] == "nanometer"  # carried over
        assert out.attributes["data_series"] == "raw_data"
        assert out.attributes["solution_a_series"] == "sol_A"
        assert table_quantity(out, CONCENTRATION_KEY).to("micromolar").magnitude == CONC
        assert table_quantity(out, PATHLENGTH_KEY).to("centimeter").magnitude == 1.0
        assert "delta epsilon" in out.name

    def test_series_default_to_most_corrected(self):
        names = [*NAMES, "a | plus_sol_A", "Background | sol_A_buffer_subtracted"]
        table = spectra_table(names, np.vstack([CD, CD[1:2], CD[:1]]))
        self.run_widget(table, titration_table([0.5]))
        assert self.widget.data_series == "plus_sol_A"
        assert self.widget.solution_a_series == "sol_A_buffer_subtracted"

    def test_saved_series_choice_is_kept(self):
        self.widget.data_series = "raw_data"
        names = [*NAMES, "a | plus_sol_A"]
        table = spectra_table(names, np.vstack([CD, CD[1:2]]))
        self.run_widget(table, titration_table([0.5, 1.0]))
        assert self.widget.data_series == "raw_data"

    def test_parameters_scale_result(self):
        self.run_widget()
        self.widget.pathlength_cm = 2.0
        self.widget.mean_residue_molecular_weight = 226.0
        self.widget.solution_a_molecular_weight = 2.0
        self.widget.commit.now()
        out = self.get_output(self.widget.Outputs.data, widget=self.widget)
        np.testing.assert_allclose(out.X[3:], expected(CD, pathlength=2.0, mrw=226.0, mw_a=2.0))

    def test_degree_rows(self):
        table = spectra_table(NAMES, CD * 1e-3, spectrum_unit="degree")
        out = self.run_widget(table)
        np.testing.assert_allclose(out.X[3:], expected(CD))
        assert spectrum_units(out)[:3] == ["degree"] * 3

    def test_concentration_unit_is_respected(self):
        titration = titration_table([0.5, 1.0], concentration=0.015)
        titration.domain["working_concentration_a"].attributes["unit"] = "millimolar"
        out = self.run_widget(titration=titration)
        np.testing.assert_allclose(out.X[3:], expected(CD))

    def test_no_input(self):
        self.send_signal(self.widget.Inputs.data, None, widget=self.widget)
        assert self.get_output(self.widget.Outputs.data, widget=self.widget) is None

    def test_without_titration_warns(self):
        self.send_signal(self.widget.Inputs.data, spectra_table(NAMES, CD), widget=self.widget)
        assert self.widget.Warning.no_titration.is_shown()
        assert self.get_output(self.widget.Outputs.data, widget=self.widget) is None

    def test_titration_without_concentration_warns(self):
        titration = Table("iris")
        out = self.run_widget(titration=titration)
        assert self.widget.Warning.missing_concentration.is_shown()
        assert out is None

    def test_concentration_without_unit_warns(self):
        titration = titration_table([0.5, 1.0])
        del titration.domain["working_concentration_a"].attributes["unit"]
        assert self.run_widget(titration=titration) is None
        assert self.widget.Warning.missing_concentration.is_shown()

    def test_invalid_parameter(self):
        self.run_widget()
        self.widget.pathlength_cm = 0.0
        self.widget.commit.now()
        assert self.widget.Error.invalid_parameter.is_shown()
        assert self.get_output(self.widget.Outputs.data, widget=self.widget) is None
        self.widget.pathlength_cm = 1.0
        self.widget.commit.now()
        assert not self.widget.Error.active

    def test_stored_invalid_parameters_are_reset(self):
        widget = self.create_widget(
            OWDeltaEpsilon,
            stored_settings=dict(pathlength_cm=-3, mean_residue_molecular_weight="x",
                                 solution_a_molecular_weight=0),
        )
        assert widget.pathlength_cm == 1.0
        assert widget.mean_residue_molecular_weight == 113.0
        assert widget.solution_a_molecular_weight == 1.0

    def test_foreign_table(self):
        assert self.run_widget(Table("iris")) is None
        assert self.widget.Error.missing_wavelength.is_shown()

    def test_names_without_stage(self):
        table = spectra_table(["x", "y", "z"], CD)
        assert self.run_widget(table) is None
        assert self.widget.Error.no_series_names.is_shown()

    def test_unknown_data_series(self):
        self.run_widget()
        self.widget.data_series = "nope"
        self.widget.commit.now()
        assert self.widget.Error.no_data_series.is_shown()

    def test_no_solution_a(self):
        self.run_widget()
        self.widget.solution_a_series = "nope"
        self.widget.commit.now()
        assert self.widget.Error.no_solution_a_series.is_shown()

    def test_missing_spectrum_unit(self):
        table = spectra_table(NAMES, CD, spectrum_unit=None)
        assert self.run_widget(table) is None
        assert "spectrum_unit" in str(self.widget.Error.invalid_parameter)

    def test_incompatible_spectrum_unit(self):
        table = spectra_table(NAMES, CD, spectrum_unit="nanometer")
        assert self.run_widget(table) is None
        assert self.widget.Error.invalid_parameter.is_shown()

    def test_concentration_label(self):
        self.run_widget()
        assert "15" in self.widget.concentration_label.text()
