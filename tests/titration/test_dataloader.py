import numpy as np
import pandas as pd
import pytest
from helpers import WAVELENGTHS, write_csv
from Orange.widgets.tests import base as orange_tests

from orangecontrib.orangeCD.widgets.titration.owcddataloader import (
    OWCDDataLoader,
    average_cd_series,
    dataframe_to_orange_table,
    file_parser,
    natural_key,
    parse_data,
    parse_remarks,
)
from orangecontrib.orangeCD.widgets.titration.owcddataloader import (
    test_empty_line as _empty_line,  # not a test
)
from orangecontrib.orangeCD.widgets.titration.utils import (
    spectrum_labels,
    spectrum_names,
    spectrum_units,
    wavelengths,
)


class TestParsers:
    def test_empty_line(self):
        assert _empty_line(",,,\n") is None
        assert _empty_line("a,b,\n") == "a, b, "

    def test_parse_remarks(self):
        remarks = parse_remarks(
            ["#Sample: lysozyme, extra\n", ",,\n", "free text,x\n", "#a:b:c\n"]
        )
        assert remarks["Sample"] == "lysozyme"
        assert remarks["Header"] == "free text, x"
        assert "a:b:c" in remarks

    def test_natural_key_sorts_numerically(self):
        names = ["/d/10 equiv.csv", "/d/2 equiv.csv", "/d/1.0 equiv.csv"]
        assert sorted(names, key=natural_key) == [
            "/d/1.0 equiv.csv",
            "/d/2 equiv.csv",
            "/d/10 equiv.csv",
        ]

    def test_natural_key_reads_decimals_as_one_number(self):
        names = ["0.3 equiv.csv", "1.05 equiv.csv", "0.25 equiv.csv", "1.5 equiv.csv"]
        assert sorted(names, key=natural_key) == [
            "0.25 equiv.csv",
            "0.3 equiv.csv",
            "1.05 equiv.csv",
            "1.5 equiv.csv",
        ]

    def test_file_parser_and_average(self, tmp_path):
        path = tmp_path / "s.csv"
        write_csv(path, [1.0, 2.0, 3.0, 4.0])
        parsed = file_parser(str(path))
        frame = parsed["Data"]["CircularDichroism"]
        assert list(frame.columns) == ["0", "1", "Average"]
        series = average_cd_series(str(path))
        assert series.name == "s.csv"
        np.testing.assert_allclose(series.index, WAVELENGTHS)
        np.testing.assert_allclose(series.values, [1, 2, 3, 4])

    def test_average_of_replicates(self, tmp_path):
        lines = [
            "Data:",
            "CircularDichroism,",
            "Wavelength,a,b",
            "nm,mdeg,mdeg",
            "260,1,3",
            "250,2,4",
            "",
            "HT,",
            "Wavelength,x",
            "nm,V",
            "260,1",
            "",
            "",
        ]
        path = tmp_path / "r.csv"
        path.write_text("\n".join(lines))
        np.testing.assert_allclose(average_cd_series(str(path)).values, [2, 3])

    def test_missing_section(self, tmp_path):
        path = tmp_path / "bad.csv"
        path.write_text("Data:\nHT,\nWavelength,x\nnm,V\n260,1\n\nHT2,\nWavelength,y\n")
        with pytest.raises(ValueError, match="CircularDichroism"):
            average_cd_series(str(path))

    def test_invalid_section(self):
        lines = ["T,", "Wavelength", "nm", "260", "", "U,", "Wavelength", ""]
        with pytest.raises(ValueError):
            parse_data(lines + ["", ""])


class TestDataframeToTable:
    def make(self):
        frame = pd.DataFrame(
            {
                ("buffer", "raw_data", "buffer"): [1.0, 2.0, 3.0],
                ("a.csv", "raw_data", "titration"): [4.0, 5.0, 6.0],
            },
            index=[260.0, 250.0, 240.0],
        )
        return dataframe_to_orange_table(frame, "angstrom", "degree")

    def test_rows_are_spectra_and_axis_ascending(self):
        table = self.make()
        assert spectrum_names(table) == ["buffer | raw_data", "a.csv | raw_data"]
        np.testing.assert_array_equal(wavelengths(table), [240, 250, 260])
        np.testing.assert_allclose(table.X, [[3, 2, 1], [6, 5, 4]])

    def test_units(self):
        table = self.make()
        assert table.attributes["wavelength_unit"] == "angstrom"
        assert spectrum_units(table) == ["degree", "degree"]


class TestLoaderWidget(orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWCDDataLoader)

    def output(self):
        return self.get_output(self.widget.Outputs.data, widget=self.widget)

    def load(self, count=2, **references):
        self.widget.data_files = [
            write_csv(self.tmp_path / f"{i}.csv", [-i] * 4) for i in range(1, count + 1)
        ]
        for setting, value in references.items():
            setattr(
                self.widget, setting, write_csv(self.tmp_path / f"{setting}.csv", value)
            )
        self.widget._refresh_file_list()
        self.widget._load_and_plot_raw_data()
        self.widget.commit()

    def test_move_selected_reorders_output(self):
        self.load(4)
        self.widget.selected_data_files = [1]
        self.widget._move_selected(-1)
        assert self.widget.data_file_names == ["2.csv", "1.csv", "3.csv", "4.csv"]
        assert list(self.widget.selected_data_files) == [0]
        assert spectrum_names(self.output())[:2] == [
            "2.csv | raw_data",
            "1.csv | raw_data",
        ]

    def test_move_selected_block_down(self):
        self.load(4)
        self.widget.selected_data_files = [1, 2]
        self.widget._move_selected(1)
        assert self.widget.data_file_names == ["1.csv", "4.csv", "2.csv", "3.csv"]
        assert list(self.widget.selected_data_files) == [2, 3]

    def test_move_selected_stops_at_the_ends(self):
        self.load(3)
        self.widget.selected_data_files = [0]
        self.widget._move_selected(-1)
        assert self.widget.data_file_names == ["1.csv", "2.csv", "3.csv"]
        self.widget.selected_data_files = [2]
        self.widget._move_selected(1)
        assert self.widget.data_file_names == ["1.csv", "2.csv", "3.csv"]

    def test_all_files_used_whatever_is_selected(self):
        self.load(3)
        assert list(self.widget.selected_data_files) == []
        assert len(spectrum_names(self.output())) == 3
        self.widget.selected_data_files = [1]
        self.widget.commit()
        assert len(spectrum_names(self.output())) == 3

    def test_empty(self):
        assert self.output() is None
        assert not self.widget.Error.active

    def test_data_files_only(self):
        self.load()
        assert spectrum_names(self.output()) == ["1.csv | raw_data", "2.csv | raw_data"]
        assert self.widget.data_file_names == ["1.csv", "2.csv"]

    def test_reference_order_and_roles(self):
        self.load(
            1,
            solution_a_file=[1] * 4,
            solution_b_file=[2] * 4,
            buffer_file=[3] * 4,
        )
        table = self.output()
        assert spectrum_names(table) == [
            "sol_A | raw_data",
            "sol_B | raw_data",
            "buffer | raw_data",
            "1.csv | raw_data",
        ]
        assert [(l.sample, l.stage, l.role) for l in spectrum_labels(table)] == [
            ("sol_A", "raw_data", "sol_A"),
            ("sol_B", "raw_data", "sol_B"),
            ("buffer", "raw_data", "buffer"),
            ("1.csv", "raw_data", "titration"),
        ]
        assert [v.name for v in table.domain.metas] == [
            "Spectrum",
            "Sample",
            "Stage",
            "Role",
        ]

    def test_references_only(self):
        self.widget.buffer_file = write_csv(self.tmp_path / "b.csv", [1] * 4)
        self.widget.commit()
        assert spectrum_names(self.output()) == ["buffer | raw_data"]

    def test_preview_has_one_curve_per_file(self):
        self.load(3)
        assert len(self.widget.plot_widget.listDataItems()) == 3

    def test_clear(self):
        self.load()
        self.widget._clear_data_files()
        assert self.widget.data_files == []
        assert self.output() is None
        assert self.widget.plot_widget.listDataItems() == []

    def test_unit_change_updates_output_and_labels(self):
        self.load()
        self.widget.wavelength_unit = "angstrom"
        self.widget.cd_unit = "degree"
        self.widget._units_changed()
        table = self.output()
        assert table.attributes["wavelength_unit"] == "angstrom"
        assert set(spectrum_units(table)) == {"degree"}
        assert "Å" in self.widget.plot_widget.getAxis("bottom").labelUnits

    def test_mismatched_wavelength_axes_use_common_wavelengths(self):
        self.widget.data_files = [
            write_csv(self.tmp_path / "a.csv", [1, 2, 3, 4]),
            write_csv(
                self.tmp_path / "b.csv",
                [5, 6, 7, 8],
                wavelengths=[250.0, 240.0, 230.0, 220.0],
            ),
        ]
        self.widget.commit()
        assert not self.widget.Error.active
        assert self.widget.Warning.cropped_wavelengths.is_shown()
        table = self.output()
        np.testing.assert_array_equal(wavelengths(table), [230, 240, 250])
        np.testing.assert_allclose(table.X, [[4, 3, 2], [7, 6, 5]])
        # The warning goes away once the axes agree.
        self.widget.data_files = self.widget.data_files[:1]
        self.widget.commit()
        assert not self.widget.Warning.cropped_wavelengths.is_shown()

    def test_no_common_wavelengths_is_an_error(self):
        self.widget.data_files = [
            write_csv(self.tmp_path / "a.csv", [1] * 4),
            write_csv(
                self.tmp_path / "b.csv",
                [1] * 4,
                wavelengths=[161.0, 151.0, 141.0, 131.0],
            ),
        ]
        self.widget.commit()
        assert self.widget.Error.load_failed.is_shown()
        assert self.output() is None

    def test_duplicate_file_names(self):
        (self.tmp_path / "x").mkdir()
        self.widget.data_files = [
            write_csv(self.tmp_path / "s.csv", [1] * 4),
            write_csv(self.tmp_path / "x" / "s.csv", [2] * 4),
        ]
        self.widget.commit()
        assert "Duplicate" in str(self.widget.Error.load_failed)
        assert self.output() is None

    def test_missing_file_clears_after_fix(self):
        self.widget.solution_a_file = str(self.tmp_path / "missing.csv")
        self.widget.commit()
        assert self.widget.Error.load_failed.is_shown()
        self.widget.solution_a_file = write_csv(self.tmp_path / "a.csv", [1] * 4)
        self.widget.commit()
        assert not self.widget.Error.load_failed.is_shown()
        assert self.output() is not None

    def test_settings_survive_reload(self):
        self.load()
        settings = self.widget.settingsHandler.pack_data(self.widget)
        widget = self.create_widget(OWCDDataLoader, stored_settings=settings)
        assert widget.data_files == self.widget.data_files
        assert spectrum_names(self.get_output(widget.Outputs.data, widget=widget)) == [
            "1.csv | raw_data",
            "2.csv | raw_data",
        ]
