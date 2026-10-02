"""Widget tests for the loader -> preprocess -> delta epsilon -> binding chain."""

import numpy as np
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.orangeCD.widgets.titration.owbindingdata import OWBindingData
from orangecontrib.orangeCD.widgets.titration.owbindingplot import OWBindingPlot
from orangecontrib.orangeCD.widgets.titration.owcddataloader import OWCDDataLoader
from orangecontrib.orangeCD.widgets.titration.owcdspectraplot import OWCDSpectraPlot
from orangecontrib.orangeCD.widgets.titration.owdeltaepsilon import OWDeltaEpsilon
from orangecontrib.orangeCD.widgets.titration.utils import (
    CONCENTRATION_KEY,
    DELTA_EPSILON_UNIT,
    spectrum_names,
    spectrum_units,
    table_quantity,
    unit_string,
    wavelengths,
)

from helpers import WAVELENGTHS, titration_table, write_csv


class LoaderFixture:
    def load(self, tmp_path):
        widget = self.create_widget(OWCDDataLoader)
        widget.solution_a_file = write_csv(tmp_path / "a.csv", [-5, -5, -5, -5])
        widget.data_files = [
            write_csv(tmp_path / f"{i}.csv", [-5 - 5 * i] * 4) for i in (1, 2, 3)
        ]
        widget._refresh_file_list()
        widget._load_and_plot_raw_data()
        widget.commit()
        return self.get_output(widget.Outputs.data, widget=widget)


class TestLoader(LoaderFixture, WidgetTest):
    def test_output_layout(self, tmp_path=None):
        import tempfile, pathlib
        with tempfile.TemporaryDirectory() as d:
            table = self.load(pathlib.Path(d))
        self.assertEqual(
            spectrum_names(table),
            ["Background | sol_A", "1.csv | raw_data", "2.csv | raw_data",
             "3.csv | raw_data"],
        )
        np.testing.assert_array_equal(wavelengths(table), sorted(WAVELENGTHS))
        self.assertEqual(table.attributes["wavelength_unit"], "nanometer")
        self.assertEqual(spectrum_units(table), ["millidegree"] * 4)

    def test_no_files_sends_nothing(self):
        widget = self.create_widget(OWCDDataLoader)
        self.assertIsNone(self.get_output(widget.Outputs.data, widget=widget))

    def test_bad_file_sets_error(self):
        widget = self.create_widget(OWCDDataLoader)
        widget.solution_a_file = "/nonexistent.csv"
        widget.commit()
        self.assertTrue(widget.Error.load_failed.is_shown())


class TestPipeline(LoaderFixture, WidgetTest):
    def run_pipeline(self):
        import tempfile, pathlib
        with tempfile.TemporaryDirectory() as d:
            loaded = self.load(pathlib.Path(d))
        # Stand-in for Quasar preprocessing: a transform keeps metas/attributes.
        from orangecontrib.spectroscopy.preprocess import Cut
        loaded = Cut(lowlim=225, highlim=265)(loaded)

        delta = self.create_widget(OWDeltaEpsilon)
        titration = titration_table([0.5, 1.0, 1.5])
        self.send_signal(delta.Inputs.titration, titration, widget=delta)
        self.send_signal(delta.Inputs.data, loaded, widget=delta)
        delta_out = self.get_output(delta.Outputs.data, widget=delta)

        binding = self.create_widget(OWBindingData)
        self.send_signal(binding.Inputs.titration, titration, widget=binding)
        self.send_signal(binding.Inputs.spectra, delta_out, widget=binding)
        return loaded, delta, delta_out, binding

    def test_delta_epsilon(self):
        _, delta, out, _ = self.run_pipeline()
        self.assertFalse(delta.Error.active)
        self.assertEqual(
            spectrum_names(out)[:4],
            ["Background | sol_A", "1.csv | raw_data", "2.csv | raw_data",
             "3.csv | raw_data"],
        )
        units = spectrum_units(out)
        self.assertEqual(units[:4], ["millidegree"] * 4)
        self.assertEqual(units[4:], [unit_string(DELTA_EPSILON_UNIT)] * 4)
        # -5 mdeg, 15 uM, 1 cm, MRW 113: -5/32980*113/15e-6
        expected = -5 / 32980 * 113 / 15e-6
        np.testing.assert_allclose(out.X[4], expected)
        self.assertAlmostEqual(
            table_quantity(out, CONCENTRATION_KEY).to("molar").magnitude, 15e-6
        )

    def test_binding_data_and_plot(self):
        _, _, delta_out, binding = self.run_pipeline()
        self.assertFalse(binding.Error.active, binding.Error.active)
        out = self.get_output(binding.Outputs.data, widget=binding)
        self.assertEqual(len(out), 4)
        self.assertEqual(out.domain["CD"].attributes["unit"], "millidegree")
        self.assertEqual(
            out.domain["Delta Epsilon"].attributes["unit"],
            unit_string(DELTA_EPSILON_UNIT),
        )
        self.assertEqual(out.domain["Conc [B]"].attributes["unit"], "molar")
        # Solution A -5 mdeg; points -10, -15, -20 -> |change| 5, 10, 15 mdeg
        np.testing.assert_allclose(out.get_column("Change in CD"), [0, 5, 10, 15])
        np.testing.assert_allclose(
            out.get_column("Delta Epsilon"),
            np.array([0, 5, 10, 15]) / 32980 / 15e-6,
        )
        plot = self.create_widget(OWBindingPlot)
        self.send_signal(plot.Inputs.data, out, widget=plot)
        self.assertEqual(plot._measurement_wavelength(), "260 nm")

    def test_spectra_plot_stage_filter(self):
        loaded, *_ = self.run_pipeline()
        widget = self.create_widget(OWCDSpectraPlot)
        self.send_signal(widget.Inputs.data, loaded, widget=widget)
        widget.selected_stage = "sol_A"
        widget._stage_changed()
        self.assertEqual(widget.spectra_names, ["Background | sol_A"])
        self.assertEqual(len(widget.plot_item.listDataItems()), 1)
        self.assertFalse(widget.Error.active)

    def test_foreign_table_shows_error_not_crash(self):
        for widget_class, signal in (
            (OWDeltaEpsilon, "data"), (OWBindingData, "spectra"),
            (OWCDSpectraPlot, "data"),
        ):
            widget = self.create_widget(widget_class)
            self.send_signal(
                getattr(widget.Inputs, signal), Table("iris"), widget=widget
            )
            self.assertTrue(widget.Error.missing_wavelength.is_shown())
