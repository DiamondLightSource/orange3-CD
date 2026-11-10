"""Pairing spectra with Titration Table rows by titration point number."""

import tempfile

import numpy as np
import pytest
from helpers import titration_table
from Orange.data import ContinuousVariable, Domain, Table
from Orange.widgets.tests import base as orange_tests
from test_correction import DILUTION, RATIO, CorrectionFixture, titration

from orangecontrib.orangeCD.widgets.titration.owbindingdata import OWBindingData
from orangecontrib.orangeCD.widgets.titration.owcddatacorrection import (
    OWCDDataCorrection,
)
from orangecontrib.orangeCD.widgets.titration.owcddataloader import OWCDDataLoader
from orangecontrib.orangeCD.widgets.titration.owdeltaepsilon import OWDeltaEpsilon
from orangecontrib.orangeCD.widgets.titration.utils import (
    POINT_META,
    build_spectra_table,
    spectrum_labels,
    spectrum_names,
    spectrum_points,
    titration_rows,
)


class TestHelpers:
    def test_titration_rows_are_zero_based(self):
        rows = titration_rows(np.array([3.0, 1.0, 2.0]), 3, "x")
        assert rows.tolist() == [2, 0, 1]

    @pytest.mark.parametrize(
        "points", [[1.0, np.nan], [0.0, 1.0], [1.0, 4.0], [1.0, 1.0]]
    )
    def test_invalid_points(self, points):
        with pytest.raises(ValueError):
            titration_rows(np.array(points), 3, "x")

    def test_points_round_trip_through_build_spectra_table(self):
        from helpers import spectra_table

        template = spectra_table(["a | s", "b | s"], np.zeros((2, 4)))
        assert spectrum_points(template) is None
        out = build_spectra_table(
            template, np.zeros((2, 4)), ["a | s", "b | s"], ["degree"] * 2,
            points=[np.nan, 2.0],
        )
        np.testing.assert_array_equal(spectrum_points(out), [np.nan, 2.0])
        assert out.domain[POINT_META].is_continuous


class PointsFixture(CorrectionFixture):
    last_loader = None

    def create_widget(self, *args, **kwargs):
        widget = super().create_widget(*args, **kwargs)
        if isinstance(widget, OWCDDataLoader):
            self.last_loader = widget
        return widget

    def load_with_titration(self, directory, titration_table_):
        """Load the files, then connect the Titration Table to the Loader."""
        self.load(directory)
        widget = self.last_loader
        self.send_signal(widget.Inputs.titration, titration_table_, widget=widget)
        return widget, self.get_output(widget.Outputs.data, widget=widget)

    @staticmethod
    def shuffled(table, order):
        out = table[order]
        out.attributes = dict(table.attributes)
        return out


class TestLoader(PointsFixture, orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWCDDataLoader)

    def test_files_are_numbered_by_list_position(self):
        with tempfile.TemporaryDirectory() as d:
            _, table = self.load_with_titration(d, titration(DILUTION, RATIO))
        names = spectrum_names(table)
        points = spectrum_points(table)
        assert points is not None
        by_name = dict(zip(names, points))
        assert [by_name[f"{i}.csv | raw_data"] for i in (1, 2, 3)] == [1, 2, 3]
        for reference in ("buffer", "sol_A", "sol_B"):
            assert np.isnan(by_name[f"{reference} | raw_data"])

    def test_list_shows_the_number_and_ratio(self):
        domain = Domain([ContinuousVariable("ratio")])
        with_ratio = Table.from_numpy(domain, RATIO[:, None])
        with tempfile.TemporaryDirectory() as d:
            widget, _ = self.load_with_titration(d, with_ratio)
        assert widget.data_file_names[0] == "1. 1.csv  (ratio 0.5)"
        assert widget.data_file_names[2] == "3. 3.csv  (ratio 1.5)"

    def test_list_shows_only_the_number_without_a_ratio_column(self):
        with tempfile.TemporaryDirectory() as d:
            widget, _ = self.load_with_titration(d, titration(DILUTION, RATIO))
        assert widget.data_file_names[0] == "1. 1.csv"

    def test_moving_a_file_renumbers_it(self):
        with tempfile.TemporaryDirectory() as d:
            widget, _ = self.load_with_titration(d, titration(DILUTION, RATIO))
            widget.selected_data_files = [2]
            widget._move_selected(-1)
            table = self.get_output(widget.Outputs.data, widget=widget)
        by_name = dict(zip(spectrum_names(table), spectrum_points(table)))
        assert by_name["3.csv | raw_data"] == 2
        assert by_name["2.csv | raw_data"] == 3

    def test_without_titration_table_there_are_no_points(self):
        with tempfile.TemporaryDirectory() as d:
            table = self.load(d)
        assert spectrum_points(table) is None

    def test_row_count_mismatch_warns_and_does_not_number(self):
        with tempfile.TemporaryDirectory() as d:
            widget, table = self.load_with_titration(
                d, titration(DILUTION[:2], RATIO[:2])
            )
        assert widget.Warning.point_count_mismatch.is_shown()
        assert spectrum_points(table) is None
        assert widget.data_file_names[0] == "1.csv"


class TestCorrectionPairing(PointsFixture, orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWCDDataCorrection)

    def corrected(self, order=None, tit=None):
        tit = titration(DILUTION, RATIO) if tit is None else tit
        with tempfile.TemporaryDirectory() as d:
            _, table = self.load_with_titration(d, tit)
        if order is not None:
            table = self.shuffled(table, order)
        return self.correct(table, tit)

    def test_row_order_does_not_matter(self):
        _, reference = self.corrected()
        with tempfile.TemporaryDirectory() as d:
            _, loaded = self.load_with_titration(d, titration(DILUTION, RATIO))
        widget, out = self.corrected(order=list(range(len(loaded)))[::-1])
        assert not widget.Error.active
        expected = dict(zip(spectrum_names(reference), reference.X))
        for name, row in zip(spectrum_names(out), out.X):
            np.testing.assert_allclose(row, expected[name], err_msg=name)

    def test_a_subset_of_points_uses_their_own_rows(self):
        _, reference = self.corrected()
        expected = dict(zip(spectrum_names(reference), reference.X))
        with tempfile.TemporaryDirectory() as d:
            _, table = self.load_with_titration(d, titration(DILUTION, RATIO))
        names = spectrum_names(table)
        keep = [i for i, n in enumerate(names) if n != "2.csv | raw_data"]
        widget, out = self.correct(
            self.shuffled(table, keep), titration(DILUTION, RATIO)
        )
        assert not widget.Error.active
        for name, row in zip(spectrum_names(out), out.X):
            if name.endswith("sol_A_added"):
                np.testing.assert_allclose(row, expected[name], err_msg=name)
        assert "2.csv | sol_A_added" not in spectrum_names(out)

    def test_points_are_carried_to_the_output(self):
        _, out = self.corrected()
        by_name = dict(zip(spectrum_names(out), spectrum_points(out)))
        assert by_name["2.csv | sol_A_added"] == 2
        assert np.isnan(by_name["sol_A | buffer_subtracted"])

    def test_point_outside_the_table_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            _, table = self.load_with_titration(d, titration(DILUTION, RATIO))
        widget, out = self.correct(table, titration(DILUTION[:2], RATIO[:2]))
        assert out is None
        assert widget.Error.invalid_titration.is_shown()


class TestDownstream(PointsFixture, orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWBindingData)

    def chain(self, reverse=False):
        tit = titration(DILUTION, RATIO)
        with tempfile.TemporaryDirectory() as d:
            _, table = self.load_with_titration(d, tit)
        correction = self.create_widget(OWCDDataCorrection)
        _, corrected = self.correct_with(correction, table, tit)
        if reverse:
            corrected = self.shuffled(corrected, list(range(len(corrected)))[::-1])
        delta = self.create_widget(OWDeltaEpsilon)
        delta.pathlength_cm = 1.0
        delta.solution_a_molecular_weight = 1.0
        self.send_signal(delta.Inputs.titration, tit, widget=delta)
        self.send_signal(delta.Inputs.data, corrected, widget=delta)
        delta_out = self.get_output(delta.Outputs.data, widget=delta)
        binding = self.create_widget(OWBindingData)
        binding.wavelength = 250.0
        self.send_signal(binding.Inputs.titration, tit, widget=binding)
        self.send_signal(binding.Inputs.spectra, delta_out, widget=binding)
        return (
            delta_out,
            binding,
            self.get_output(binding.Outputs.data, widget=binding),
        )

    def correct_with(self, widget, table, tit):
        self.send_signal(widget.Inputs.titration, tit, widget=widget)
        self.send_signal(widget.Inputs.data, table, widget=widget)
        return widget, self.get_output(widget.Outputs.data, widget=widget)

    def test_binding_data_is_in_point_order_whatever_the_row_order(self):
        _, _, reference = self.chain()
        assert reference is not None
        _, binding, shuffled = self.chain(reverse=True)
        assert not binding.Error.active
        np.testing.assert_allclose(shuffled.X, reference.X)
        assert list(shuffled.metas[:, 0]) == list(reference.metas[:, 0])
        np.testing.assert_allclose(shuffled.get_column("Titration point"), [0, *RATIO])

    def test_delta_epsilon_carries_points(self):
        delta_out, _, _ = self.chain()
        by_name = dict(zip(spectrum_names(delta_out), spectrum_points(delta_out)))
        assert by_name["3.csv | sol_A_added"] == 3
        assert by_name["3.csv | sol_A_added_delta_epsilon"] == 3

    def test_binding_data_rejects_unknown_points(self):
        _, binding, _ = self.chain()
        self.send_signal(
            binding.Inputs.titration, titration_table(RATIO[:2]), widget=binding
        )
        assert binding.Error.mismatched_points.is_shown()


class TestSeparatorInFileName(PointsFixture, orange_tests.WidgetTest):
    """A file name containing ' | ' must not be misread as sample | stage."""

    def setUp(self):
        self.widget = self.create_widget(OWCDDataCorrection)

    def test_pipeline_keeps_the_sample_and_roles(self):
        from pathlib import Path

        from helpers import write_csv
        from test_correction import RAW, REFERENCES, instrument_order

        with tempfile.TemporaryDirectory() as d:
            directory = Path(d)
            loader = self.create_widget(OWCDDataLoader)
            loader.buffer_file = write_csv(
                directory / "buffer.csv", instrument_order(REFERENCES["buffer"])
            )
            loader.solution_a_file = write_csv(
                directory / "a.csv", instrument_order(REFERENCES["sol_A"])
            )
            loader.data_files = [
                write_csv(directory / f"x | {i}.csv", instrument_order(raw))
                for i, raw in enumerate(RAW, 1)
            ]
            loader._refresh_file_list()
            loader._load_and_plot_raw_data()
            loader.commit()
            loaded = self.get_output(loader.Outputs.data, widget=loader)

        correction = self.create_widget(OWCDDataCorrection)
        self.send_signal(
            correction.Inputs.titration, titration(DILUTION, RATIO), widget=correction
        )
        self.send_signal(correction.Inputs.data, loaded, widget=correction)
        out = self.get_output(correction.Outputs.data, widget=correction)
        assert not correction.Error.active, correction.Error.active
        labels = spectrum_labels(out)
        added = [l for l in labels if l.stage == "sol_A_added"]
        assert [l.sample for l in added] == ["x | 1.csv", "x | 2.csv", "x | 3.csv"]
        assert {l.role for l in added} == {"titration"}
        assert ("sol_A", "buffer_subtracted", "sol_A") in [
            (l.sample, l.stage, l.role) for l in labels
        ]
