from Orange.widgets.tests.base import WidgetTest

from orangecontrib.orangeCD.widgets.titration.owtitration import OWTitrationCalculator


class TestTitrationVolumeWarning(WidgetTest):
    def make(self, mode, ratios, **settings):
        widget = self.create_widget(OWTitrationCalculator)
        widget.mode = mode
        widget.ratios = ratios
        for name, value in settings.items():
            setattr(widget, name, value)
        widget.calculate()
        return widget

    def test_increasing_within_limit_has_no_warning(self):
        widget = self.make("increasing", "0.14, 0.28, 0.42")
        self.assertFalse(widget.Warning.volume_exceeded.is_shown())
        self.assertFalse(widget.Error.active)

    def test_increasing_exceeding_limit_warns(self):
        widget = self.make(
            "increasing", "1.6, 2.71, 3.14, 6.282, 10, 20, 40, 80",
            stock_b_concentrations_values="200",
        )
        self.assertTrue(widget.Warning.volume_exceeded.is_shown())
        self.assertFalse(widget.Error.active)
        self.assertIsNotNone(self.get_output(widget.Outputs.data, widget=widget))

    def test_warning_clears_when_recalculated_within_limit(self):
        widget = self.make(
            "increasing", "1.6, 3.14, 10, 20, 40, 80", stock_b_concentrations_values="200",
        )
        self.assertTrue(widget.Warning.volume_exceeded.is_shown())
        widget.ratios = "0.14, 0.28"
        widget.stock_b_concentrations_values = "2000, 4000, 8000"
        widget.calculate()
        self.assertFalse(widget.Warning.volume_exceeded.is_shown())

    def test_fixed_mode_never_warns(self):
        widget = self.make("fixed", "1.6, 3.14", stock_b_concentrations_values="200")
        self.assertFalse(widget.Warning.volume_exceeded.is_shown())


class TestTitrationOutputColumns(WidgetTest):
    def output(self, mode, ratios):
        widget = self.create_widget(OWTitrationCalculator)
        widget.mode = mode
        widget.ratios = ratios
        widget.calculate()
        return self.get_output(widget.Outputs.data, widget=widget)

    def test_increasing_columns_and_attributes(self):
        table = self.output("increasing", "0.14, 0.28, 0.42")
        self.assertEqual(
            [v.name for v in table.domain.attributes],
            ["ratio", "stock_b", "volume_added_this_step", "total_stock_b_volume",
             "total_cell_volume", "dilution_factor", "normalised_molar_ratio",
             "working_concentration_a"],
        )
        self.assertEqual(len(table.domain.metas), 0)
        self.assertEqual(table.domain["working_concentration_a"].attributes["unit"], "micromolar")
        self.assertEqual(table.attributes["mode"], "increasing")
        self.assertTrue(table.attributes["within_limit"])
        for key in ("volume_solution_a", "volume_buffer", "max_volume_allowed", "max_volume_added"):
            self.assertIn("microliter", table.attributes[key])

    def test_fixed_columns_and_attributes(self):
        table = self.output("fixed", "0.5, 1.0")
        self.assertEqual(
            [v.name for v in table.domain.attributes],
            ["ratio", "stock_b", "volume_stock_b", "baseline_volume",
             "concentration_b", "normalised_molar_ratio", "working_concentration_a"],
        )
        self.assertEqual(table.attributes["mode"], "fixed")
        self.assertNotIn("within_limit", table.attributes)
