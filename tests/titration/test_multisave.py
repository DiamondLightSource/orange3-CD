from unittest import mock

import numpy as np
import pandas as pd
import pytest
from Orange.data import (
    ContinuousVariable,
    DiscreteVariable,
    Domain,
    StringVariable,
    Table,
)
from Orange.widgets.tests import base as orange_tests

from orangecontrib.orangeCD.widgets.titration import owmultisave
from orangecontrib.orangeCD.widgets.titration.owmultisave import (
    OWMultiSave,
    make_unique_sheet_names,
    table_to_dataframe,
)


def small_table(name="t", offset=0.0):
    domain = Domain(
        [ContinuousVariable("a"), ContinuousVariable("b")],
        DiscreteVariable("cls", values=("x", "y")),
        metas=[StringVariable("label")],
    )
    table = Table.from_numpy(
        domain,
        np.array([[1.0, 2.0], [3.0, 4.0]]) + offset,
        np.array([0, 1]),
        metas=np.array([["r1"], ["r2"]], dtype=object),
    )
    table.name = name
    return table


class TestSheetNames:
    def test_plain_names_are_kept(self):
        assert make_unique_sheet_names(["a", "b"]) == ["a", "b"]

    def test_duplicates_ignoring_case(self):
        assert make_unique_sheet_names(["Data", "data", "DATA"]) == [
            "Data",
            "data (1)",
            "DATA (2)",
        ]

    def test_invalid_characters_are_replaced(self):
        assert make_unique_sheet_names(["a/b\\c*d?e[f]g:h"]) == ["a_b_c_d_e_f_g_h"]

    def test_long_names_are_truncated_and_stay_unique(self):
        names = make_unique_sheet_names(["x" * 40, "x" * 40])
        assert all(len(name) <= 31 for name in names)
        assert len(set(names)) == 2

    @pytest.mark.parametrize("empty", ["", None, "   ", "///"])
    def test_empty_names(self, empty):
        (name,) = make_unique_sheet_names([empty])
        assert name.startswith(("Sheet", "___"))
        assert name


class TestTableToDataframe:
    def test_columns_and_values(self):
        frame = table_to_dataframe(small_table())
        assert list(frame.columns) == ["a", "b", "cls", "label"]
        assert frame["a"].tolist() == [1.0, 3.0]
        assert frame["label"].tolist() == ["r1", "r2"]
        assert frame["cls"].tolist() == [0.0, 1.0]

    def test_attributes_only(self):
        table = Table.from_numpy(
            Domain([ContinuousVariable("a")]), np.array([[1.0], [2.0]])
        )
        assert list(table_to_dataframe(table).columns) == ["a"]


class TestHeaderUnits:
    def test_variable_units_are_added(self):
        conc = ContinuousVariable("Conc [B]")
        conc.attributes["unit"] = "molar"
        delta = ContinuousVariable("Delta Epsilon")
        delta.attributes["unit"] = "liter / mole / centimeter"
        point = ContinuousVariable("Titration point")  # no unit
        ratio = ContinuousVariable("dilution_factor")
        ratio.attributes["unit"] = "dimensionless"  # nothing worth showing
        table = Table.from_numpy(
            Domain([point, conc, delta, ratio]), np.ones((2, 4))
        )
        assert list(table_to_dataframe(table).columns) == [
            "Titration point",
            "Conc [B] (M)",
            "Delta Epsilon (M\u207b\u00b9 cm\u207b\u00b9)",
            "dilution_factor",
        ]

    def test_wavelength_columns_take_the_table_unit(self):
        table = Table.from_numpy(
            Domain([ContinuousVariable("230.0"), ContinuousVariable("240.0")]),
            np.ones((1, 2)),
        )
        table.attributes["wavelength_unit"] = "nanometer"
        assert list(table_to_dataframe(table).columns) == ["230.0 (nm)", "240.0 (nm)"]

    def test_non_numeric_names_do_not_take_the_wavelength_unit(self):
        table = Table.from_numpy(
            Domain([ContinuousVariable("a")]), np.ones((1, 1))
        )
        table.attributes["wavelength_unit"] = "nanometer"
        assert list(table_to_dataframe(table).columns) == ["a"]

    def test_unit_pint_does_not_know_is_shown_as_written(self):
        variable = ContinuousVariable("x")
        variable.attributes["unit"] = "widgets per fortnight"
        table = Table.from_numpy(Domain([variable]), np.ones((1, 1)))
        assert list(table_to_dataframe(table).columns) == ["x (widgets per fortnight)"]


class TestEmptyTable:
    def test_table_without_columns(self):
        table = Table.from_domain(Domain([]), 3)
        frame = table_to_dataframe(table)
        assert len(frame) == 3
        assert len(frame.columns) == 0


class TestWidget(orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWMultiSave)

    def send(self, table, index):
        self.send_signal(self.widget.Inputs.data, table, index, widget=self.widget)

    def save(self, filename):
        target = "orangecontrib.orangeCD.widgets.titration.owmultisave.QFileDialog.getSaveFileName"
        with mock.patch(target, return_value=(filename, "")):
            self.widget.save_workbook()

    def test_inputs_become_sheets_named_after_tables(self):
        self.send(small_table("first"), 0)
        self.send(small_table("second"), 1)
        assert [item["sheet_name"] for item in self.widget.model.items] == [
            "first",
            "second",
        ]

    def test_default_sheet_name_when_table_unnamed(self):
        table = small_table()
        table.name = ""
        self.send(table, 0)
        assert self.widget.model.items[0]["sheet_name"] == "Sheet1"

    def test_update_and_remove(self):
        self.send(small_table("first"), 0)
        self.send(small_table("changed", 10), 0)
        assert self.widget.model.items[0]["table"].X[0, 0] == 11.0
        self.send_signal(
            self.widget.Inputs.data,
            self.widget.Inputs.data.closing_sentinel,
            0,
            widget=self.widget,
        )
        assert self.widget.model.items == []

    def test_sheet_name_is_editable(self):
        self.send(small_table("first"), 0)
        index = self.widget.model.index(0, 1)
        assert self.widget.model.setData(index, "renamed", owmultisave.Qt.EditRole)
        assert self.widget.model.data(index, owmultisave.Qt.DisplayRole) == "renamed"
        assert not self.widget.model.setData(
            self.widget.model.index(0, 0), "x", owmultisave.Qt.EditRole
        )

    def test_save_workbook(self):
        path = self.tmp_path / "book"
        self.send(small_table("Data"), 0)
        self.send(small_table("data", 10), 1)
        self.save(str(path))
        sheets = pd.read_excel(str(path) + ".xlsx", sheet_name=None)
        assert list(sheets) == ["Data", "data (1)"]
        assert sheets["Data"]["a"].tolist() == [1.0, 3.0]
        assert sheets["data (1)"]["a"].tolist() == [11.0, 13.0]
        assert sheets["Data"]["label"].tolist() == ["r1", "r2"]

    def test_table_attributes_go_to_a_metadata_sheet(self):
        first = small_table("Data")
        first.attributes.update({"solution_a_concentration": "15.0 micromolar"})
        second = small_table("Other")
        second.attributes.update({"within_limit": True, "region": [1.0, 2.0]})
        self.send(first, 0)
        self.send(second, 1)
        self.send(small_table("Plain"), 2)
        self.save(str(self.tmp_path / "book"))
        sheets = pd.read_excel(str(self.tmp_path / "book.xlsx"), sheet_name=None)
        assert list(sheets) == ["Data", "Other", "Plain", "Metadata"]
        meta = sheets["Metadata"]
        assert list(meta.columns) == ["Sheet", "Attribute", "Value"]
        assert meta.values.tolist() == [
            ["Data", "solution_a_concentration", "15.0 micromolar"],
            ["Other", "within_limit", "True"],
            ["Other", "region", "[1.0, 2.0]"],
        ]

    def test_no_metadata_sheet_without_attributes(self):
        self.send(small_table("Data"), 0)
        self.save(str(self.tmp_path / "book"))
        sheets = pd.read_excel(str(self.tmp_path / "book.xlsx"), sheet_name=None)
        assert list(sheets) == ["Data"]

    def test_metadata_sheet_name_does_not_clash(self):
        table = small_table("Metadata")
        table.attributes["k"] = "v"
        self.send(table, 0)
        self.save(str(self.tmp_path / "book"))
        sheets = pd.read_excel(str(self.tmp_path / "book.xlsx"), sheet_name=None)
        assert list(sheets) == ["Metadata", "Metadata (1)"]

    def test_nothing_to_save_shows_error(self):
        self.save(str(self.tmp_path / "book"))
        assert self.widget.Error.save_failed.is_shown()
        assert list(self.tmp_path.iterdir()) == []
        self.send(small_table("Data"), 0)
        self.save(str(self.tmp_path / "book"))
        assert not self.widget.Error.save_failed.is_shown()

    def test_unwritable_destination_shows_error(self):
        self.send(small_table("Data"), 0)
        self.save(str(self.tmp_path / "missing_folder" / "book"))
        assert self.widget.Error.save_failed.is_shown()

    def test_extension_is_not_duplicated(self):
        self.send(small_table("Data"), 0)
        self.save(str(self.tmp_path / "book.xlsx"))
        assert [p.name for p in self.tmp_path.iterdir()] == ["book.xlsx"]

    def test_cancelled_dialog_writes_nothing(self):
        self.send(small_table(), 0)
        self.save("")
        assert list(self.tmp_path.iterdir()) == []
