"""Tests for the CD Data Correction widget and its use in the pipeline."""

import pathlib
import tempfile

import numpy as np
from Orange.data import ContinuousVariable, Domain, Table
from Orange.widgets.tests import base as orange_tests

from orangecontrib.orangeCD.widgets.titration.owbindingdata import OWBindingData
from orangecontrib.orangeCD.widgets.titration.owcddatacorrection import (
    OWCDDataCorrection,
    correct_spectra,
)
from orangecontrib.orangeCD.widgets.titration.owcddataloader import OWCDDataLoader
from orangecontrib.orangeCD.widgets.titration.owdeltaepsilon import OWDeltaEpsilon
from orangecontrib.orangeCD.widgets.titration.utils import (
    spectrum_names,
    spectrum_units,
)

from helpers import WAVELENGTHS, write_csv

ASC = sorted(WAVELENGTHS)
SHAPE = np.arange(4.0)  # wavelength-dependent part, in ascending order


def instrument_order(values):
    """Spectra are written in descending wavelength order."""
    return list(np.asarray(values)[::-1])


def titration(dilution, ratio, concentration=15.0):
    domain = Domain([
        ContinuousVariable("dilution_factor"),
        ContinuousVariable("normalised_molar_ratio"),
        ContinuousVariable("working_concentration_a"),
    ])
    domain.attributes[2].attributes["unit"] = "micromolar"
    return Table.from_numpy(
        domain,
        np.column_stack((dilution, ratio, np.full(len(ratio), concentration))),
    )


REFERENCES = {
    "buffer": 1 + 0.1 * SHAPE,
    "sol_A": 6 + 0.2 * SHAPE,
    "sol_B": 3 + 0.3 * SHAPE,
}
RAW = [10 + 2 * i + 0.5 * SHAPE for i in range(3)]
DILUTION = np.array([1.0, 1.05, 1.1])
RATIO = np.array([0.5, 1.0, 1.5])


class CorrectionFixture:
    def load(self, directory, with_sol_b=True):
        widget = self.create_widget(OWCDDataLoader)
        d = pathlib.Path(directory)
        widget.buffer_file = write_csv(d / "buffer.csv", instrument_order(REFERENCES["buffer"]))
        widget.solution_a_file = write_csv(d / "a.csv", instrument_order(REFERENCES["sol_A"]))
        if with_sol_b:
            widget.solution_b_file = write_csv(d / "b.csv", instrument_order(REFERENCES["sol_B"]))
        widget.data_files = [
            write_csv(d / f"{i}.csv", instrument_order(raw)) for i, raw in enumerate(RAW, 1)
        ]
        widget._refresh_file_list()
        widget._load_and_plot_raw_data()
        widget.commit()
        return self.get_output(widget.Outputs.data, widget=widget)

    def correct(self, table, titration_table):
        widget = self.widget
        self.send_signal(widget.Inputs.titration, titration_table, widget=widget)
        self.send_signal(widget.Inputs.data, table, widget=widget)
        return widget, self.get_output(widget.Outputs.data, widget=widget)


class TestCorrection(CorrectionFixture, orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWCDDataCorrection)

    def test_matches_reference_calculation(self):
        with tempfile.TemporaryDirectory() as d:
            loaded = self.load(d)
        widget, out = self.correct(loaded, titration(DILUTION, RATIO))
        self.assertFalse(widget.Error.active)
        names = spectrum_names(out)
        by_name = dict(zip(names, out.X))

        buffer, sol_a, sol_b = (REFERENCES[k] for k in ("buffer", "sol_A", "sol_B"))
        a_bs, b_bs = sol_a - buffer, sol_b - buffer
        np.testing.assert_allclose(by_name["Background | sol_A_buffer_subtracted"], a_bs)
        np.testing.assert_allclose(by_name["Background | sol_B_buffer_subtracted"], b_bs)
        for i, raw in enumerate(RAW):
            name = f"{i + 1}.csv"
            buffer_sub = (raw - buffer) * DILUTION[i]
            sol_a_sub = buffer_sub - a_bs
            frac = sol_a_sub - b_bs * RATIO[i]
            np.testing.assert_allclose(by_name[f"{name} | buffer_subtraction"], buffer_sub)
            np.testing.assert_allclose(by_name[f"{name} | sol_A_subtraction"], sol_a_sub)
            np.testing.assert_allclose(by_name[f"{name} | subtract_frac_sol_B"], frac)
            np.testing.assert_allclose(by_name[f"{name} | plus_sol_A"], frac + a_bs)
            # inputs are retained untouched
            np.testing.assert_allclose(by_name[f"{name} | raw_data"], raw)
        self.assertEqual(set(spectrum_units(out)), {"millidegree"})
        self.assertEqual(out.attributes["wavelength_unit"], "nanometer")

    def test_unit_conversion_of_references(self):
        with tempfile.TemporaryDirectory() as d:
            loaded = self.load(d)
        # Rewrite the buffer row in degrees: 1e-3 * mdeg values.
        names = spectrum_names(loaded)
        row = names.index("Background | buffer")
        x = loaded.X.copy()
        x[row] *= 1e-3
        from orangecontrib.orangeCD.widgets.titration.utils import build_spectra_table
        units = ["millidegree"] * len(names)
        units[row] = "degree"
        table = build_spectra_table(loaded, x, names, units)
        _, out = self.correct(table, titration(DILUTION, RATIO))
        by_name = dict(zip(spectrum_names(out), out.X))
        np.testing.assert_allclose(
            by_name["Background | sol_A_buffer_subtracted"],
            REFERENCES["sol_A"] - REFERENCES["buffer"],
        )

    def test_without_sol_b_warns_and_skips(self):
        with tempfile.TemporaryDirectory() as d:
            loaded = self.load(d, with_sol_b=False)
        widget, out = self.correct(loaded, titration(DILUTION, RATIO))
        self.assertTrue(widget.Warning.no_sol_b.is_shown())
        by_name = dict(zip(spectrum_names(out), out.X))
        np.testing.assert_allclose(
            by_name["1.csv | subtract_frac_sol_B"], by_name["1.csv | sol_A_subtraction"]
        )

    def test_errors(self):
        with tempfile.TemporaryDirectory() as d:
            loaded = self.load(d)
        widget = self.widget
        self.send_signal(widget.Inputs.data, loaded, widget=widget)
        self.assertTrue(widget.Warning.no_titration.is_shown())
        self.assertIsNone(self.get_output(widget.Outputs.data, widget=widget))
        self.send_signal(widget.Inputs.titration, titration(DILUTION[:2], RATIO[:2]), widget=widget)
        self.assertTrue(widget.Error.invalid_titration.is_shown())
        self.send_signal(widget.Inputs.titration, Table("iris"), widget=widget)
        self.assertTrue(widget.Error.invalid_titration.is_shown())
        self.send_signal(widget.Inputs.data, Table("iris"), widget=widget)
        self.assertTrue(widget.Error.missing_spectrum_meta.is_shown())

    def test_missing_dilution_factor_assumes_one_and_warns(self):
        with tempfile.TemporaryDirectory() as d:
            loaded = self.load(d)
        no_dilution = Table.from_numpy(
            Domain([ContinuousVariable("normalised_molar_ratio")]),
            RATIO[:, None],
        )
        widget, out = self.correct(loaded, no_dilution)
        self.assertTrue(widget.Warning.no_dilution_factor.is_shown())
        self.assertFalse(widget.Error.active)
        self.assertIsNotNone(out)
        by_name = dict(zip(spectrum_names(out), out.X))
        np.testing.assert_allclose(
            by_name["2.csv | buffer_subtraction"], RAW[1] - REFERENCES["buffer"]
        )
        # a titration table with the column clears the warning
        self.send_signal(widget.Inputs.titration, titration(DILUTION, RATIO), widget=widget)
        self.assertFalse(widget.Warning.no_dilution_factor.is_shown())

    def test_missing_background(self):
        with tempfile.TemporaryDirectory() as d:
            widget = self.create_widget(OWCDDataLoader)
            widget.data_files = [write_csv(pathlib.Path(d) / "1.csv", [1, 2, 3, 4])]
            widget._refresh_file_list()
            widget._load_and_plot_raw_data()
            widget.commit()
            loaded = self.get_output(widget.Outputs.data, widget=widget)
        corr, _ = self.correct(loaded, titration([1.0], [0.5]))
        self.assertTrue(corr.Error.missing_background.is_shown())

    def test_full_chain_defaults_to_corrected_series(self):
        with tempfile.TemporaryDirectory() as d:
            loaded = self.load(d)
        _, corrected = self.correct(loaded, titration(DILUTION, RATIO))
        tit = titration(DILUTION, RATIO)
        delta = self.create_widget(OWDeltaEpsilon)
        self.send_signal(delta.Inputs.titration, tit, widget=delta)
        self.send_signal(delta.Inputs.data, corrected, widget=delta)
        self.assertEqual(delta.data_series, "plus_sol_A")
        self.assertEqual(delta.solution_a_series, "sol_A_buffer_subtracted")
        delta_out = self.get_output(delta.Outputs.data, widget=delta)
        self.assertFalse(delta.Error.active)

        binding = self.create_widget(OWBindingData)
        self.send_signal(binding.Inputs.titration, tit, widget=binding)
        self.send_signal(binding.Inputs.spectra, delta_out, widget=binding)
        self.assertEqual(binding.data_series, "plus_sol_A")
        self.assertFalse(binding.Error.active, binding.Error.active)
        out = self.get_output(binding.Outputs.data, widget=binding)
        self.assertEqual(len(out), 4)


def test_correct_spectra_without_sol_b():
    result = correct_spectra(
        np.array([[3.0, 4.0]]), np.array([1.0, 1.0]), np.array([2.0, 2.0]),
        None, np.array([2.0]), np.array([0.5]),
    )
    np.testing.assert_allclose(result["buffer_subtraction"], [[4.0, 6.0]])
    np.testing.assert_allclose(result["plus_sol_A"], [[4.0, 6.0]])
