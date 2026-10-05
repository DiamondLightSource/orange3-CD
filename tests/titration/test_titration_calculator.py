import pytest
from Orange.widgets.tests import base as orange_tests
from pint import Quantity

from orangecontrib.orangeCD.widgets.titration import Q_
from orangecontrib.orangeCD.widgets.titration.owtitration import (
    CDTitrationCalculator,
    OWTitrationCalculator,
    TitrationMode,
    TitrationPoint,
)

UL = "microliter"
UM = "micromolar"


def make_calculator(**overrides):
    arguments = dict(
        starting_cell_volume=Q_(500.0, UL),
        stock_con_a=Q_(468.0, UM),
        working_con_a=Q_(19.659, UM),
        stock_b_concentrations=[Q_(2000, UM), Q_(4000, UM), Q_(8000, UM)],
        stock_b_molar_equiv=2.75,
    )
    arguments.update(overrides)
    return CDTitrationCalculator(**arguments)


def magnitude(quantity, unit=UL):
    return float(quantity.to(unit).magnitude)


class TestConstruction:
    @pytest.mark.parametrize(
        "override",
        [
            dict(starting_cell_volume=Q_(0, UL)),
            dict(stock_con_a=Q_(0, UM)),
            dict(working_con_a=Q_(-1, UM)),
            dict(stock_b_concentrations=[]),
            dict(stock_b_concentrations=[Q_(1, UM), Q_(0, UM)]),
        ],
    )
    def test_invalid(self, override):
        with pytest.raises(ValueError):
            make_calculator(**override)

    def test_volume_of_solution_a(self):
        assert magnitude(make_calculator().volume_solution_a) == 21.0

    def test_volume_follows_units(self):
        calculator = make_calculator(starting_cell_volume=Q_(0.5, "milliliter"))
        assert magnitude(calculator.volume_solution_a) == pytest.approx(21.0, abs=0.1)


class TestChooseStock:
    def test_prefers_least_concentrated_in_range(self):
        stock, volume = make_calculator().choose_stock(1.6)
        assert stock == 1
        assert magnitude(volume) == pytest.approx(1.6 * 19.659 * 500 / 2000)

    def test_skips_stock_that_is_too_dilute(self):
        stock, _ = make_calculator().choose_stock(5.0)  # stock 1 -> 24.6 uL
        assert stock == 2

    def test_nothing_in_range_picks_closest(self):
        assert make_calculator().choose_stock(0.14)[0] == 1  # all below the minimum
        assert make_calculator().choose_stock(20.0)[0] == 3  # all above the maximum

    def test_custom_range(self):
        stock, _ = make_calculator().choose_stock(
            0.14, target_min=Q_(0.1, UL), target_max=Q_(1.0, UL)
        )
        assert stock == 1

    @pytest.mark.parametrize(
        "ratio, low, high",
        [(0, 2, 20), (-1, 2, 20), (1, -1, 20), (1, 5, 5), (1, 10, 5)],
    )
    def test_invalid(self, ratio, low, high):
        with pytest.raises(ValueError):
            make_calculator().choose_stock(
                ratio, target_min=Q_(low, UL), target_max=Q_(high, UL)
            )


class TestCreatePoints:
    def test_fixed_uses_full_ratio(self):
        points = make_calculator().create_points([1.6, 3.2], mode=TitrationMode.FIXED)
        assert [p.ratio for p in points] == [1.6, 3.2]
        assert magnitude(points[1].predicted_volume) == pytest.approx(3.2 * 19.659 * 500 / 2000)

    def test_increasing_uses_increment(self):
        points = make_calculator().create_points([1.6, 2.6], mode=TitrationMode.INCREASING)
        assert magnitude(points[1].predicted_volume) == pytest.approx(1.0 * 19.659 * 500 / 2000)

    @pytest.mark.parametrize("ratios", [[], [0.5, 0], [-1]])
    def test_invalid_ratios(self, ratios):
        with pytest.raises(ValueError):
            make_calculator().create_points(ratios, mode=TitrationMode.FIXED)

    def test_increasing_must_increase(self):
        with pytest.raises(ValueError, match="increasing"):
            make_calculator().create_points([1, 1], mode=TitrationMode.INCREASING)
        # the same ratios are fine for fixed volumes
        make_calculator().create_points([2, 1], mode=TitrationMode.FIXED)

    def test_unknown_mode(self):
        with pytest.raises(ValueError):
            make_calculator().create_points([1], mode="sideways")


class TestCalculateFixed:
    def setup_method(self):
        calculator = make_calculator()
        points = [TitrationPoint(ratio=1.0, stock_b=1), TitrationPoint(ratio=2.0, stock_b=2)]
        self.result = calculator.calculate(TitrationMode.FIXED, points)

    def test_rows(self):
        first, second = self.result.rows
        assert magnitude(first.volume_stock_b) == 4.9
        assert magnitude(first.baseline_volume) == pytest.approx(474.1)
        assert magnitude(first.concentration_b, UM) == 19.7
        assert first.normalised_molar_ratio == 0.364
        assert magnitude(second.volume_stock_b) == 4.9  # stock 2 is twice as strong
        assert magnitude(second.baseline_volume) == pytest.approx(474.1)

    def test_cell_volume_is_constant(self):
        for row in self.result.rows:
            total = row.volume_stock_b + row.baseline_volume + self.result.volume_solution_a
            assert magnitude(total) == pytest.approx(500.0)

    def test_dataframe_has_only_per_point_columns(self):
        frame, units = self.result.result_to_dataframe()
        assert list(frame.columns) == [
            "ratio", "stock_b", "volume_stock_b", "baseline_volume",
            "concentration_b", "normalised_molar_ratio", "working_concentration_a",
        ]
        assert units["baseline_volume"] == UL
        assert units["working_concentration_a"] == UM
        assert "ratio" not in units
        assert not any(isinstance(v, Quantity) for v in frame.to_numpy().ravel())

    def test_attributes(self):
        attributes = self.result.result_attributes()
        assert attributes["mode"] == "fixed"
        assert set(attributes) == {"mode", "volume_solution_a"}

    def test_negative_baseline(self):
        calculator = make_calculator(
            starting_cell_volume=Q_(30, UL),
            stock_b_concentrations=[Q_(100, UM)],
        )
        with pytest.raises(ValueError, match="negative"):
            calculator.calculate(TitrationMode.FIXED, [TitrationPoint(ratio=10.0, stock_b=1)])


class TestCalculateIncreasing:
    def setup_method(self):
        calculator = make_calculator()
        points = calculator.create_points([0.14, 0.28], mode=TitrationMode.INCREASING)
        self.result = calculator.calculate(TitrationMode.INCREASING, points)

    def test_running_totals(self):
        first, second = self.result.rows
        assert magnitude(first.volume_added_this_step) == 0.7
        assert magnitude(second.total_stock_b_volume) == 1.4
        assert magnitude(second.total_cell_volume) == pytest.approx(501.4)
        assert first.dilution_factor == 1.001
        assert second.dilution_factor == 1.003

    def test_limits(self):
        assert magnitude(self.result.volume_buffer) == 479.0
        assert magnitude(self.result.max_volume_allowed) == 75.0
        assert magnitude(self.result.max_volume_added) == 1.4
        assert self.result.within_limit

    def test_attributes(self):
        attributes = self.result.result_attributes()
        assert attributes["mode"] == "increasing"
        assert attributes["within_limit"] is True
        assert attributes["max_volume_allowed"].endswith("microliter")

    def test_exceeding_limit(self):
        calculator = make_calculator(stock_b_concentrations=[Q_(100, UM)])
        points = [TitrationPoint(ratio=r, stock_b=1) for r in (5.0, 10.0, 15.0)]
        result = calculator.calculate(TitrationMode.INCREASING, points)
        assert not result.within_limit
        assert result.result_attributes()["within_limit"] is False


class TestCalculateValidation:
    @pytest.mark.parametrize(
        "points",
        [[], [TitrationPoint(0.0, 1)], [TitrationPoint(1.0, 0)], [TitrationPoint(1.0, 4)]],
    )
    def test_invalid_points(self, points):
        with pytest.raises(ValueError):
            make_calculator().calculate(TitrationMode.FIXED, points)

    def test_unknown_mode(self):
        with pytest.raises(ValueError):
            make_calculator().calculate("sideways", [TitrationPoint(1.0, 1)])


class TestCalculatorWidget(orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWTitrationCalculator)

    def calculate(self, **settings):
        for name, value in settings.items():
            setattr(self.widget, name, value)
        self.widget.calculate()
        return self.get_output(self.widget.Outputs.data, widget=self.widget)

    def test_default_run(self):
        table = self.calculate()
        assert table is not None
        assert not self.widget.Error.active
        assert self.widget.table_model.rowCount() == 6

    def test_invalid_input_shows_error_and_clears_output(self):
        for settings in (
            dict(ratios="abc"),
            dict(ratios=""),
            dict(ratios="2, 1", mode="increasing"),
            dict(stock_b_concentrations_values="x"),
            dict(starting_cell_volume_value="-5"),
            dict(stock_con_a_value="zero"),
        ):
            with self.subTest(**settings):
                self.setUp()
                table = self.calculate(**settings)
                assert table is None
                assert self.widget.Error.invalid_input.is_shown()
                assert self.widget.table_model.rowCount() == 0

    def test_error_clears_after_valid_input(self):
        self.calculate(ratios="abc")
        assert self.widget.Error.invalid_input.is_shown()
        assert self.calculate(ratios="1, 2") is not None
        assert not self.widget.Error.invalid_input.is_shown()

    def test_unit_change_updates_labels_and_output(self):
        self.calculate(ratios="1.0", mode="fixed")
        self.widget._conc_unit = "nanomolar"
        self.widget._volume_unit = "milliliter"
        self.widget._on_unit_changed()
        symbols = [label.text() for label, *_ in self.widget._unit_labels]
        assert any("(ml)" in text for text in symbols)
        assert any("nM" in text for text in symbols)
        table = self.calculate()
        assert table.domain["baseline_volume"].attributes["unit"] == "milliliter"
        assert table.domain["working_concentration_a"].attributes["unit"] == "nanomolar"

    def test_summary_mentions_limit_in_increasing_mode(self):
        self.calculate(mode="increasing", ratios="0.14, 0.28")
        assert "within 15% limit: yes" in self.widget.summary_label.text()
        self.calculate(mode="fixed", ratios="1.0")
        assert "limit" not in self.widget.summary_label.text()
