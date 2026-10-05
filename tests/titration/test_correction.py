"""Tests for the CD Data Correction widget and its use in the pipeline."""

import pathlib
import tempfile

import numpy as np
import pytest
from helpers import WAVELENGTHS, write_csv
from Orange.data import ContinuousVariable, Domain, Table
from Orange.widgets.tests import base as orange_tests

from orangecontrib.orangeCD.widgets.titration.owbindingdata import OWBindingData
from orangecontrib.orangeCD.widgets.titration.owcddatacorrection import (
    OWCDDataCorrection,
    auto_baseline,
    correct_spectra,
    find_flat_region,
)
from orangecontrib.orangeCD.widgets.titration.owcddataloader import OWCDDataLoader
from orangecontrib.orangeCD.widgets.titration.owdeltaepsilon import OWDeltaEpsilon
from orangecontrib.orangeCD.widgets.titration.utils import (
    spectrum_names,
    spectrum_units,
)

ASC = sorted(WAVELENGTHS)
SHAPE = np.arange(4.0)  # wavelength-dependent part, in ascending order


def instrument_order(values):
    """Spectra are written in descending wavelength order."""
    return list(np.asarray(values)[::-1])


def titration(dilution, ratio, concentration=15.0):
    domain = Domain(
        [
            ContinuousVariable("dilution_factor"),
            ContinuousVariable("normalised_molar_ratio"),
            ContinuousVariable("working_concentration_a"),
        ]
    )
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
        widget.buffer_file = write_csv(
            d / "buffer.csv", instrument_order(REFERENCES["buffer"])
        )
        widget.solution_a_file = write_csv(
            d / "a.csv", instrument_order(REFERENCES["sol_A"])
        )
        if with_sol_b:
            widget.solution_b_file = write_csv(
                d / "b.csv", instrument_order(REFERENCES["sol_B"])
            )
        widget.data_files = [
            write_csv(d / f"{i}.csv", instrument_order(raw))
            for i, raw in enumerate(RAW, 1)
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
        np.testing.assert_allclose(
            by_name["Background | sol_A_buffer_subtracted"], a_bs
        )
        np.testing.assert_allclose(
            by_name["Background | sol_B_buffer_subtracted"], b_bs
        )
        for i, raw in enumerate(RAW):
            name = f"{i + 1}.csv"
            buffer_sub = (raw - buffer) * DILUTION[i]
            sol_a_sub = buffer_sub - a_bs
            frac = sol_a_sub - b_bs * RATIO[i]
            np.testing.assert_allclose(
                by_name[f"{name} | buffer_subtraction"], buffer_sub
            )
            np.testing.assert_allclose(
                by_name[f"{name} | sol_A_subtraction"], sol_a_sub
            )
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
        self.send_signal(
            widget.Inputs.titration, titration(DILUTION[:2], RATIO[:2]), widget=widget
        )
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
        self.send_signal(
            widget.Inputs.titration, titration(DILUTION, RATIO), widget=widget
        )
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
        np.array([[3.0, 4.0]]),
        np.array([1.0, 1.0]),
        np.array([2.0, 2.0]),
        None,
        np.array([2.0]),
        np.array([0.5]),
    )
    np.testing.assert_allclose(result["buffer_subtraction"], [[4.0, 6.0]])
    np.testing.assert_allclose(result["plus_sol_A"], [[4.0, 6.0]])


WL = np.arange(190.0, 261.0)  # 71 points, signal below 230 nm, flat above 235 nm
LAZY_NAMES = [
    "Background | buffer",
    "Background | sol_A",
    "Background | sol_B",
    "a | raw_data",
    "b | raw_data",
    "c | raw_data",
]


def synthetic_spectra(
    offsets=(0.3, 1.0, 0.5, 0.7, 0.9, 1.2),
    amplitudes=(0.1, 6.0, 3.0, 2.0, 4.0, 5.0),
    noise=0.002,
):
    rng = np.random.default_rng(0)
    signal = np.exp(-(((WL - 210) / 8) ** 2)) * np.sin((WL - 190) / 6)
    return np.array(
        [
            offset + amplitude * signal + rng.normal(0, noise, WL.size)
            for offset, amplitude in zip(offsets, amplitudes)
        ]
    )


def lazy_table(X=None, **kwargs):
    from helpers import spectra_table

    return spectra_table(
        LAZY_NAMES, synthetic_spectra() if X is None else X, wavelengths=WL, **kwargs
    )


class TestFlatRegion:
    def test_finds_the_flat_stretch(self):
        start, stop = find_flat_region(synthetic_spectra(), 7)
        assert stop - start == 7
        assert WL[start] >= 235  # beyond the signal band

    def test_ignores_a_constant_offset_and_scale(self):
        base = synthetic_spectra(offsets=(0,) * 6, amplitudes=(1,) * 6)
        scaled = base * np.array([1, 50, 0.1, 3, 7, 2])[:, None] + 100
        assert find_flat_region(base, 7) == find_flat_region(scaled, 7)

    def test_signal_at_the_other_end(self):
        X = synthetic_spectra()[:, ::-1]
        _, stop = find_flat_region(X, 7)
        assert stop <= 36  # the flat part is now at the start of the axis

    def test_window_is_clamped(self):
        X = synthetic_spectra()
        assert find_flat_region(X, 10_000) == (0, X.shape[1])
        start, stop = find_flat_region(X, 0)
        assert stop - start == 1

    def test_all_flat_spectra(self):
        assert find_flat_region(np.ones((3, 20)), 5) == (0, 5)

    def test_nan_values_do_not_break_the_search(self):
        X = synthetic_spectra()
        X[1, 5] = np.nan
        start, _ = find_flat_region(X, 7)
        assert WL[start] >= 235


class TestAutoBaseline:
    def test_subtracts_the_mean_of_the_flat_region(self):
        X = synthetic_spectra()
        corrected, offsets, (start, stop) = auto_baseline(X, 10)
        np.testing.assert_allclose(offsets, X[:, start:stop].mean(axis=1))
        np.testing.assert_allclose(corrected, X - offsets[:, None])
        np.testing.assert_allclose(corrected[:, start:stop].mean(axis=1), 0, atol=1e-12)

    def test_recovers_the_known_offsets(self):
        _, offsets, _ = auto_baseline(synthetic_spectra(), 10)
        np.testing.assert_allclose(offsets, [0.3, 1.0, 0.5, 0.7, 0.9, 1.2], atol=0.01)

    def test_window_width_follows_percentage(self):
        X = synthetic_spectra()
        _, _, (start, stop) = auto_baseline(X, 30)
        assert stop - start == round(X.shape[1] * 0.30)

    def test_minimum_window(self):
        _, _, (start, stop) = auto_baseline(synthetic_spectra(), 0)
        assert stop - start == 3

    def test_input_is_not_modified(self):
        X = synthetic_spectra()
        before = X.copy()
        auto_baseline(X)
        np.testing.assert_array_equal(X, before)


class TestPreprocessingWarnings(CorrectionFixture, orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWCDDataCorrection)
        self.tit = titration([1.0, 1.0, 1.0], RATIO)

    def send(self, table, tit=None):
        self.send_signal(
            self.widget.Inputs.titration,
            self.tit if tit is None else tit,
            widget=self.widget,
        )
        self.send_signal(self.widget.Inputs.data, table, widget=self.widget)

    def preprocessed(self, table):
        preprocess = pytest.importorskip("orangecontrib.spectroscopy.preprocess")
        return preprocess.LinearBaseline()(table)

    def test_raw_input_warns_but_still_emits(self):
        self.send(lazy_table())
        assert self.widget.Warning.no_preprocessing.is_shown()
        assert not self.widget.Warning.double_baseline.is_shown()
        assert self.get_output(self.widget.Outputs.data, widget=self.widget) is not None
        assert not self.widget.Error.active

    def test_preprocessed_input_does_not_warn(self):
        self.send(self.preprocessed(lazy_table()))
        assert not self.widget.Warning.no_preprocessing.is_shown()
        assert self.get_output(self.widget.Outputs.data, widget=self.widget) is not None

    def test_warning_clears_when_preprocessed_data_arrives(self):
        table = lazy_table()
        self.send(table)
        assert self.widget.Warning.no_preprocessing.is_shown()
        self.send(self.preprocessed(table))
        assert not self.widget.Warning.no_preprocessing.is_shown()

    def test_warning_is_independent_of_the_titration_table(self):
        self.send_signal(self.widget.Inputs.data, lazy_table(), widget=self.widget)
        assert self.widget.Warning.no_preprocessing.is_shown()
        assert self.widget.Warning.no_titration.is_shown()

    def test_no_data_means_no_warning(self):
        self.send(lazy_table())
        self.send_signal(self.widget.Inputs.data, None, widget=self.widget)
        assert not self.widget.Warning.no_preprocessing.is_shown()
        assert self.get_output(self.widget.Outputs.data, widget=self.widget) is None

    def test_foreign_table_reports_the_error_not_a_crash(self):
        self.send(Table("iris"))
        assert self.widget.Error.missing_spectrum_meta.is_shown()


class TestLazyProcess(CorrectionFixture, orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWCDDataCorrection)
        self.tit = titration([1.0, 1.0, 1.0], RATIO)

    def send(self, table):
        self.send_signal(self.widget.Inputs.titration, self.tit, widget=self.widget)
        self.send_signal(self.widget.Inputs.data, table, widget=self.widget)

    def output(self):
        return self.get_output(self.widget.Outputs.data, widget=self.widget)

    def tick(self, checked=True):
        from AnyQt.QtWidgets import QCheckBox

        box = next(
            c
            for c in self.widget.controlArea.findChildren(QCheckBox)
            if c.text() == "Lazy process"
        )
        if box.isChecked() != checked:
            box.click()
        return box

    def rows(self, table):
        return dict(zip(spectrum_names(table), table.X))

    def test_checkbox_exists_and_is_off_by_default(self):
        box = self.tick(checked=False)
        assert not box.isChecked() and not self.widget.lazy_process

    def test_ticking_enables_the_setting_and_recalculates(self):
        self.send(lazy_table())
        before = self.rows(self.output())["a | plus_sol_A"]
        box = self.tick()
        assert box.isChecked() and self.widget.lazy_process
        after = self.rows(self.output())["a | plus_sol_A"]
        assert not np.allclose(before, after)

    def test_unticking_restores_the_plain_result(self):
        self.send(lazy_table())
        before = self.rows(self.output())["a | plus_sol_A"]
        self.tick()
        self.tick(checked=False)
        np.testing.assert_allclose(self.rows(self.output())["a | plus_sol_A"], before)

    def test_matches_manual_baseline_subtraction(self):
        from orangecontrib.orangeCD.widgets.titration.utils import build_spectra_table

        X = synthetic_spectra()
        manual_X, _, _ = auto_baseline(X, self.widget.lazy_window_percent)
        manual = build_spectra_table(
            lazy_table(X), manual_X, LAZY_NAMES, ["millidegree"] * 6
        )
        self.send(manual)
        expected = self.output().X[6:]
        self.send(lazy_table(X))
        self.tick()
        np.testing.assert_allclose(self.output().X[6:], expected, atol=1e-12)

    def test_input_rows_are_kept_unmodified(self):
        X = synthetic_spectra()
        self.send(lazy_table(X))
        self.tick()
        np.testing.assert_allclose(self.output().X[:6], X)

    def test_baseline_is_removed_from_the_corrected_spectra(self):
        # Constant offsets must not leak into the corrected result.
        flat = synthetic_spectra(amplitudes=(0.1, 6.0, 3.0, 2.0, 4.0, 5.0))
        self.send(lazy_table(flat))
        self.tick()
        out = self.rows(self.output())
        tail = WL >= 236
        for stage in (
            "buffer_subtraction",
            "sol_A_subtraction",
            "subtract_frac_sol_B",
            "plus_sol_A",
        ):
            assert np.abs(out[f"a | {stage}"][tail]).max() < 0.05, stage

    def test_region_label_and_attribute(self):
        self.send(lazy_table())
        assert self.widget.lazy_label.text() == ""
        self.tick()
        start, stop = (
            float(self.output().attributes["lazy_baseline_region"][i]) for i in (0, 1)
        )
        assert start >= 235 and stop <= 260
        assert f"{start:g}" in self.widget.lazy_label.text()
        assert f"{stop:g}" in self.widget.lazy_label.text()
        assert "nm" in self.widget.lazy_label.text()

    def test_no_region_attribute_without_lazy_process(self):
        self.send(lazy_table())
        assert "lazy_baseline_region" not in self.output().attributes

    def test_unticking_clears_the_label(self):
        self.send(lazy_table())
        self.tick()
        self.tick(checked=False)
        assert self.widget.lazy_label.text() == ""

    def test_ticking_silences_the_no_preprocessing_warning(self):
        self.send(lazy_table())
        assert self.widget.Warning.no_preprocessing.is_shown()
        self.tick()
        assert not self.widget.Warning.no_preprocessing.is_shown()
        assert not self.widget.Warning.double_baseline.is_shown()
        self.tick(checked=False)
        assert self.widget.Warning.no_preprocessing.is_shown()

    def test_double_baseline_warning_for_preprocessed_input(self):
        preprocess = pytest.importorskip("orangecontrib.spectroscopy.preprocess")
        self.send(preprocess.LinearBaseline()(lazy_table()))
        self.tick()
        assert self.widget.Warning.double_baseline.is_shown()
        assert not self.widget.Warning.no_preprocessing.is_shown()
        assert self.output() is not None  # still emitted
        self.tick(checked=False)
        assert not self.widget.Warning.double_baseline.is_shown()

    def test_window_width_changes_the_region(self):
        self.send(lazy_table())
        self.tick()
        narrow = self.output().attributes["lazy_baseline_region"]
        self.widget.lazy_window_percent = 40
        self.widget.commit()
        wide = self.output().attributes["lazy_baseline_region"]
        assert (wide[1] - wide[0]) > (narrow[1] - narrow[0])

    def test_settings_survive_reload(self):
        self.tick()
        self.widget.lazy_window_percent = 25
        settings = self.widget.settingsHandler.pack_data(self.widget)
        widget = self.create_widget(OWCDDataCorrection, stored_settings=settings)
        assert widget.lazy_process is True
        assert widget.lazy_window_percent == 25

    def test_errors_still_apply_in_lazy_mode(self):
        self.tick()
        table = lazy_table()
        # drop the buffer row: the corrections cannot run
        names = [n for n in LAZY_NAMES if n != "Background | buffer"]
        from helpers import spectra_table

        self.send(spectra_table(names, table.X[1:], wavelengths=WL))
        assert self.widget.Error.missing_background.is_shown()
        assert self.output() is None
        assert self.widget.lazy_label.text() == ""

    def test_works_with_other_wavelength_units(self):
        self.tick()
        self.send(lazy_table(wavelength_unit="angstrom"))
        assert "\u00c5" in self.widget.lazy_label.text()

    def test_no_wavelength_unit(self):
        self.tick()
        self.send(lazy_table(wavelength_unit=None))
        assert self.widget.lazy_label.text().startswith(
            "Baseline taken from the flat region"
        )
        assert self.output() is not None
