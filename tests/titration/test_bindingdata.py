import numpy as np
from helpers import spectra_table, titration_table
from Orange.data import Table
from Orange.widgets.tests import base as orange_tests

from orangecontrib.orangeCD.widgets.titration import Q_
from orangecontrib.orangeCD.widgets.titration.owbindingdata import OWBindingData
from orangecontrib.orangeCD.widgets.titration.utils import (
    CONCENTRATION_KEY,
    MEASUREMENT_WAVELENGTH_KEY,
    PATHLENGTH_KEY,
    quantity_string,
    table_quantity,
)

NAMES = ["Background | sol_A", "a | raw_data", "b | raw_data", "c | raw_data"]
# columns are 230, 240, 250, 260 nm
CD = np.array(
    [
        [0.0, 0.0, -5.0, 0.0],
        [1.0, 1.0, -10.0, 1.0],
        [2.0, 2.0, -15.0, 2.0],
        [3.0, 3.0, -20.0, 3.0],
    ]
)
RATIOS = [0.5, 1.0, 1.5]
CONC_M = 15e-6


def spectra(**kwargs):
    table = spectra_table(NAMES, CD, **kwargs)
    table.attributes[CONCENTRATION_KEY] = quantity_string(Q_(15, "micromolar"))
    table.attributes[PATHLENGTH_KEY] = quantity_string(Q_(1, "centimeter"))
    return table


class TestBindingData(orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWBindingData)
        self.widget.wavelength = 250.0

    def run_widget(self, table=None, titration=None):
        self.send_signal(
            self.widget.Inputs.titration,
            titration_table(RATIOS) if titration is None else titration,
            widget=self.widget,
        )
        self.send_signal(
            self.widget.Inputs.spectra,
            spectra() if table is None else table,
            widget=self.widget,
        )
        return self.get_output(self.widget.Outputs.data, widget=self.widget)

    def column(self, out, name):
        return out.get_column(name)

    def test_absolute_change_values(self):
        out = self.run_widget()
        assert not self.widget.Error.active
        np.testing.assert_allclose(self.column(out, "CD"), [-5, -10, -15, -20])
        np.testing.assert_allclose(self.column(out, "Change in CD"), [0, 5, 10, 15])
        delta_a = np.array([0, 5, 10, 15]) / 32980
        np.testing.assert_allclose(self.column(out, "Delta A"), delta_a)
        np.testing.assert_allclose(self.column(out, "Delta Epsilon"), delta_a / CONC_M)
        titration_point = np.array([0, *RATIOS])
        np.testing.assert_allclose(self.column(out, "Titration point"), titration_point)
        np.testing.assert_allclose(
            self.column(out, "Binding Stoichiometry"),
            titration_point / (titration_point + 1),
        )
        np.testing.assert_allclose(
            self.column(out, "Conc [B]"), titration_point * CONC_M
        )
        assert list(out.metas[:, 0]) == NAMES

    def test_signed_change(self):
        self.widget.absolute_change = False
        out = self.run_widget()
        np.testing.assert_allclose(self.column(out, "Change in CD"), [0, -5, -10, -15])
        assert np.all(self.column(out, "Delta Epsilon") <= 0)

    def test_units_on_output_columns(self):
        out = self.run_widget()
        units = {v.name: v.attributes.get("unit") for v in out.domain.attributes}
        assert units["CD"] == "millidegree"
        assert units["Change in CD"] == "millidegree"
        assert units["Conc [B]"] == "molar"
        assert units["Delta Epsilon"] == "liter / centimeter / mole"
        assert units["Delta A"] is None  # dimensionless

    def test_degree_input_gives_same_physics(self):
        table = spectra_table(NAMES, CD * 1e-3, spectrum_unit="degree")
        table.attributes.update(
            {
                key: value
                for key, value in spectra().attributes.items()
                if key in (CONCENTRATION_KEY, PATHLENGTH_KEY)
            }
        )
        out = self.run_widget(table)
        assert out.domain["CD"].attributes["unit"] == "degree"
        np.testing.assert_allclose(self.column(out, "CD"), CD[:, 2] * 1e-3)
        np.testing.assert_allclose(
            self.column(out, "Delta Epsilon"),
            np.array([0, 5, 10, 15]) / 32980 / CONC_M,
        )

    def test_mixed_row_units_are_converted(self):
        units = ["millidegree", "degree", "millidegree", "millidegree"]
        X = CD.copy()
        X[1] *= 1e-3  # "a" given in degrees
        table = spectra_table(NAMES, X, units=units, spectrum_unit=None)
        table.attributes.update(
            {
                key: value
                for key, value in spectra().attributes.items()
                if key in (CONCENTRATION_KEY, PATHLENGTH_KEY)
            }
        )
        out = self.run_widget(table)
        np.testing.assert_allclose(self.column(out, "CD"), [-5, -10, -15, -20])
        assert out.domain["CD"].attributes["unit"] == "millidegree"

    def test_metadata(self):
        out = self.run_widget()
        assert out.attributes["absolute_change"] is True
        assert out.attributes["data_series"] == "raw_data"
        assert (
            table_quantity(out, MEASUREMENT_WAVELENGTH_KEY).to("nanometer").magnitude
            == 250
        )
        assert table_quantity(out, CONCENTRATION_KEY).to("micromolar").magnitude == 15
        assert out.name == "Binding and Origin data at 250 nm"

    def test_wavelength_snaps_to_nearest_column(self):
        self.widget.wavelength = 253.0
        out = self.run_widget()
        assert self.widget.wavelength == 250.0
        np.testing.assert_allclose(self.column(out, "CD"), CD[:, 2])
        self.widget.wavelength = 100.0
        self.widget._wavelength_control_changed()
        assert self.widget.wavelength == 230.0
        out = self.get_output(self.widget.Outputs.data, widget=self.widget)
        np.testing.assert_allclose(self.column(out, "CD"), CD[:, 0])

    def test_dragging_the_line_snaps_and_commits(self):
        self.run_widget()
        self.widget.line.setValue(239.0)
        self.widget._line_changed()
        assert self.widget.wavelength == 240.0
        out = self.get_output(self.widget.Outputs.data, widget=self.widget)
        np.testing.assert_allclose(self.column(out, "CD"), CD[:, 1])

    def test_plot_shows_reference_and_series(self):
        self.run_widget()
        assert len(self.widget.plot.listDataItems()) == 4
        assert self.widget.plot.getAxis("bottom").labelUnits == "nm"
        assert self.widget.plot.getAxis("left").labelUnits == "mdeg"
        assert self.widget.wavelength_spin.suffix() == " nm"

    def test_plot_axis_without_units(self):
        self.run_widget(
            spectra_table(NAMES, CD, wavelength_unit=None, spectrum_unit=None)
        )
        assert self.widget.wavelength_spin.suffix() == ""

    def test_series_controls_and_changes(self):
        names = [*NAMES, "a | plus_sol_A", "b | plus_sol_A", "c | plus_sol_A"]
        X = np.vstack([CD, CD[1:] * 2])
        table = spectra_table(names, X)
        table.attributes.update(spectra().attributes)
        out = self.run_widget(table)
        assert self.widget.data_series == "plus_sol_A"
        np.testing.assert_allclose(self.column(out, "CD"), [-5, -20, -30, -40])
        self.widget.data_series = "raw_data"
        self.widget._series_changed()
        out = self.get_output(self.widget.Outputs.data, widget=self.widget)
        np.testing.assert_allclose(self.column(out, "CD"), [-5, -10, -15, -20])

    def test_nothing_without_both_inputs(self):
        self.send_signal(self.widget.Inputs.spectra, spectra(), widget=self.widget)
        assert self.get_output(self.widget.Outputs.data, widget=self.widget) is None
        assert not self.widget.Error.active

    def test_foreign_spectra(self):
        assert self.run_widget(Table("iris")) is None
        assert self.widget.Error.missing_wavelength.is_shown()

    def test_invalid_wavelength_axis(self):
        table = spectra_table(NAMES, CD, wavelengths=("a", "b", "c", "d"))
        assert self.run_widget(table) is None
        assert self.widget.Error.invalid_wavelength.is_shown()

    def test_missing_titration_ratios(self):
        assert self.run_widget(titration=Table("iris")) is None
        assert self.widget.Error.missing_ratios.is_shown()

    def test_mismatched_point_count(self):
        assert self.run_widget(titration=titration_table([0.5, 1.0])) is None
        assert "2 titration ratios" in str(self.widget.Error.mismatched_points)

    def test_missing_conversion_metadata(self):
        table = spectra()
        del table.attributes[CONCENTRATION_KEY]
        assert self.run_widget(table) is None
        assert self.widget.Error.missing_conversion_metadata.is_shown()

    def test_invalid_conversion_values(self):
        for value in ("0 centimeter", "nonsense", "-2 centimeter"):
            with self.subTest(value=value):
                self.setUp()
                table = spectra()
                table.attributes[PATHLENGTH_KEY] = value
                assert self.run_widget(table) is None
                assert self.widget.Error.invalid_conversion_metadata.is_shown()

    def test_missing_reference(self):
        self.run_widget()
        with self.assertWarns(UserWarning):  # Orange warns about the stale value
            self.widget.solution_a_series = "nope"
        self.widget.commit()
        assert self.widget.Error.missing_solution_a.is_shown()

    def test_missing_data_series(self):
        self.run_widget()
        with self.assertWarns(UserWarning):
            self.widget.data_series = "nope"
        self.widget.commit()
        assert self.widget.Error.missing_data.is_shown()

    def test_delta_epsilon_series_is_rejected(self):
        units = ["millidegree"] + ["liter / mole / centimeter"] * 3
        table = spectra_table(NAMES, CD, units=units, spectrum_unit=None)
        table.attributes.update(spectra().attributes)
        table.attributes.pop("spectrum_unit", None)
        out = self.run_widget(table)
        assert out is None
        assert self.widget.Error.incompatible_units.is_shown()

    def test_missing_spectrum_unit(self):
        table = spectra()
        del table.attributes["spectrum_unit"]
        assert self.run_widget(table) is None
        assert self.widget.Error.incompatible_units.is_shown()

    def test_errors_clear_on_recovery(self):
        self.run_widget(titration=Table("iris"))
        assert self.widget.Error.missing_ratios.is_shown()
        self.send_signal(
            self.widget.Inputs.titration, titration_table(RATIOS), widget=self.widget
        )
        assert not self.widget.Error.active
