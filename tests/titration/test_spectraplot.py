import numpy as np
import pytest
from AnyQt.QtGui import QColor
from helpers import spectra_table
from Orange.data import Table
from Orange.widgets.tests import base as orange_tests

from orangecontrib.orangeCD.widgets.titration.owcdspectraplot import (
    COLOUR_SCALES,
    OWCDSpectraPlot,
    colours_from_scale,
)

NAMES = ["Background | sol_A", "a | raw_data", "b | raw_data", "a | plus_sol_A"]
X = np.arange(16.0).reshape(4, 4)


class TestColours:
    def test_zero(self):
        assert colours_from_scale("Viridis", 0) == []

    def test_single_colour_is_middle_anchor(self):
        (colour,) = colours_from_scale("Classic", 1)
        assert colour == QColor(COLOUR_SCALES["Classic"][3])

    @pytest.mark.parametrize("scale", list(COLOUR_SCALES))
    def test_ends_match_anchors(self, scale):
        colours = colours_from_scale(scale, 7)
        assert len(colours) == 7
        anchors = COLOUR_SCALES[scale]
        assert colours[0].name() == QColor(anchors[0]).name()
        assert colours[-1].name() == QColor(anchors[-1]).name()

    def test_interpolates(self):
        colours = colours_from_scale("Greyscale", 3)
        assert colours[0].red() < colours[1].red() < colours[2].red()


class TestNormaliseSetting:
    choices = ("A", "B")

    def test_index_to_text(self):
        assert (
            OWCDSpectraPlot._normalise_combo_setting(1, self.choices, default="A")
            == "B"
        )

    def test_bad_index_and_value(self):
        normalise = OWCDSpectraPlot._normalise_combo_setting
        assert normalise(9, self.choices, default="A") == "A"
        assert normalise("C", self.choices, default="B") == "B"
        assert normalise("B", self.choices, default="A") == "B"


class TestWidget(orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWCDSpectraPlot)

    def send(self, table):
        self.send_signal(self.widget.Inputs.data, table, widget=self.widget)

    def curves(self):
        return self.widget.plot_item.listDataItems()

    def legend_labels(self):
        return [label.text for _, label in self.widget.legend.items]

    def test_stages_sorted_with_all_first(self):
        self.send(spectra_table(NAMES, X))
        assert self.widget.processing_stages == [
            "All spectra",
            "plus_sol_A",
            "raw_data",
            "sol_A",
        ]

    def test_all_spectra_selected_and_plotted(self):
        self.send(spectra_table(NAMES, X))
        assert self.widget.spectra_names == NAMES
        assert self.widget.selected_spectra == [0, 1, 2, 3]
        assert len(self.curves()) == 4
        assert not self.widget.Error.active

    def test_curve_data_matches_rows(self):
        self.send(spectra_table(NAMES, X))
        x, y = self.curves()[2].getData()
        np.testing.assert_array_equal(x, [230, 240, 250, 260])
        np.testing.assert_array_equal(y, X[2])

    def test_stage_filter_plots_the_right_rows(self):
        self.send(spectra_table(NAMES, X))
        self.widget.selected_stage = "raw_data"
        self.widget._stage_changed()
        assert self.widget.spectra_names == ["a | raw_data", "b | raw_data"]
        ys = [curve.getData()[1].tolist() for curve in self.curves()]
        assert ys == [X[1].tolist(), X[2].tolist()]
        assert self.legend_labels() == ["a", "b"]

    def test_partial_selection(self):
        self.send(spectra_table(NAMES, X))
        self.widget.selected_stage = "raw_data"
        self.widget._stage_changed()
        self.widget.selected_spectra = [1]
        self.widget._replot()
        (curve,) = self.curves()
        np.testing.assert_array_equal(curve.getData()[1], X[2])
        assert self.legend_labels() == ["b"]

    def test_select_all_button(self):
        self.send(spectra_table(NAMES, X))
        self.widget.selected_spectra = [0]
        self.widget._select_all_spectra()
        assert len(self.curves()) == 4

    def test_stage_is_kept_when_new_data_has_it(self):
        self.send(spectra_table(NAMES, X))
        self.widget.selected_stage = "raw_data"
        self.widget._stage_changed()
        self.send(spectra_table(NAMES, X))
        assert self.widget.selected_stage == "raw_data"
        self.send(spectra_table(["x | other"] * 4, X))
        assert self.widget.selected_stage == "All spectra"

    def test_empty_selection_plots_nothing(self):
        self.send(spectra_table(NAMES, X))
        self.widget.selected_spectra = []
        self.widget._replot()
        assert self.curves() == []
        assert not self.widget.Error.no_numeric_spectra.is_shown()

    def test_legend_toggle(self):
        self.send(spectra_table(NAMES, X))
        assert self.legend_labels() == ["Background", "a", "b", "a"]
        self.widget.show_legend = False
        self.widget._replot()
        assert self.legend_labels() == []

    def test_reverse_axis(self):
        self.send(spectra_table(NAMES, X))
        assert not self.widget.plot_item.vb.xInverted()
        self.widget.reverse_wavelength_axis = True
        self.widget._replot()
        assert self.widget.plot_item.vb.xInverted()

    def test_non_finite_points_are_dropped(self):
        data = X.copy()
        data[1, 2] = np.nan
        self.send(spectra_table(NAMES, data))
        x, y = self.curves()[1].getData()
        assert len(x) == len(y) == 3

    def test_colour_scale_changes_pen(self):
        self.send(spectra_table(NAMES, X))
        before = self.curves()[0].opts["pen"].color().name()
        self.widget.colour_scale = "Plasma"
        self.widget._replot()
        assert self.curves()[0].opts["pen"].color().name() != before

    def test_line_width(self):
        self.send(spectra_table(NAMES, X))
        self.widget.line_width = 4.0
        self.widget._replot()
        assert self.curves()[0].opts["pen"].widthF() == 4.0

    def test_axis_units(self):
        self.send(spectra_table(NAMES, X))
        assert self.widget.plot_item.getAxis("bottom").labelUnits == "nm"
        assert self.widget.plot_item.getAxis("left").labelUnits == "mdeg"

    def test_y_axis_unit_hidden_when_units_mixed(self):
        units = ["millidegree", "millidegree", "degree", "millidegree"]
        self.send(spectra_table(NAMES, X, units=units))
        assert self.widget.plot_item.getAxis("left").labelUnits == ""
        self.widget.selected_stage = "raw_data"
        self.widget._stage_changed()
        assert self.widget.plot_item.getAxis("left").labelUnits == ""
        self.widget.selected_spectra = [0]
        self.widget._replot()
        assert self.widget.plot_item.getAxis("left").labelUnits == "mdeg"

    def test_no_units(self):
        self.send(spectra_table(NAMES, X, spectrum_unit=None, wavelength_unit=None))
        assert self.widget.plot_item.getAxis("bottom").labelUnits == ""

    def test_foreign_table(self):
        self.send(Table("iris"))
        assert self.widget.Error.missing_wavelength.is_shown()
        assert self.curves() == []

    def test_invalid_wavelengths(self):
        self.send(spectra_table(NAMES, X, wavelengths=("a", "b", "c", "d")))
        assert self.widget.Error.invalid_wavelength.is_shown()

    def test_clearing_input(self):
        self.send(spectra_table(NAMES, X))
        self.send(None)
        assert self.curves() == []
        assert self.widget.spectra_names == []
        assert self.widget.processing_stages == ["All spectra"]

    def test_empty_table_reports_no_spectra(self):
        self.send(spectra_table(NAMES, X)[:0])
        assert self.widget.Error.no_numeric_spectra.is_shown()
