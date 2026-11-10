import numpy as np
import pytest
from helpers import spectra_table
from Orange.data import ContinuousVariable, Domain, Table

from orangecontrib.orangeCD.units import Q_
from orangecontrib.orangeCD.widgets.titration.utils import (
    MDEG_PER_DELTA_A,
    InvalidWavelength,
    MissingWavelength,
    ROLE_BUFFER,
    ROLE_SOL_A,
    ROLE_SOL_B,
    ROLE_TITRATION,
    STAGE_RAW,
    SpectraError,
    SpectrumLabel,
    build_spectra_table,
    has_spectroscopy_preprocessing,
    matching_spectra,
    quantity_string,
    shared_unit,
    spectrum_labels,
    spectrum_names,
    spectrum_units,
    split_series_name,
    stages,
    table_quantity,
    table_unit,
    unit_string,
    unit_symbol,
    wavelengths,
)

NAMES = ["Background | sol_A", "a | raw_data", "b | raw_data", "a | sol_A_added"]


@pytest.fixture
def table():
    return spectra_table(NAMES, np.arange(16.0).reshape(4, 4))


class TestSeriesNames:
    @pytest.mark.parametrize(
        "name, expected",
        [
            ("sample | raw_data", ("sample", "raw_data")),
            ("  x  |  stage ", ("x", "stage")),
            ("a | b | c", ("a | b", "c")),
            ("no separator", None),
            (" | stage", None),
            ("sample | ", None),
        ],
    )
    def test_split(self, name, expected):
        assert split_series_name(name) == expected


class TestLayout:
    def test_names_and_wavelengths(self, table):
        assert spectrum_names(table) == NAMES
        np.testing.assert_array_equal(wavelengths(table), [230, 240, 250, 260])

    def test_missing_meta(self):
        with pytest.raises(MissingWavelength):
            spectrum_names(Table("iris"))
        with pytest.raises(SpectraError):
            wavelengths(Table("iris"))

    def test_meta_must_be_string(self):
        domain = Domain(
            [ContinuousVariable("1.0")], metas=[ContinuousVariable("Spectrum")]
        )
        table = Table.from_numpy(domain, np.zeros((1, 1)), metas=np.zeros((1, 1)))
        with pytest.raises(MissingWavelength):
            spectrum_names(table)

    def test_invalid_wavelength(self):
        table = spectra_table(["a | x"], [[1, 2]], wavelengths=("1.0", "abc"))
        with pytest.raises(InvalidWavelength):
            wavelengths(table)

    def test_matching_spectra(self, table):
        assert matching_spectra(table) == list(enumerate(NAMES))
        # the Solution A reference is a raw_data spectrum of its own role
        assert matching_spectra(table, "raw_data") == [
            (0, "Background | sol_A"),
            (1, "a | raw_data"),
            (2, "b | raw_data"),
        ]
        assert matching_spectra(table, "raw_data", ROLE_TITRATION) == [
            (1, "a | raw_data"),
            (2, "b | raw_data"),
        ]
        assert matching_spectra(table, "raw_data", ROLE_SOL_A) == [
            (0, "Background | sol_A")
        ]
        assert matching_spectra(table, "missing") == []
        # the stage must match exactly, not just end the name
        assert matching_spectra(table, "data") == []

    def test_stages_in_order_of_appearance(self, table):
        assert stages(table) == ["raw_data", "sol_A_added"]
        assert stages(table, ROLE_SOL_A) == ["raw_data"]
        assert stages(table, ROLE_TITRATION) == ["raw_data", "sol_A_added"]


class TestLabels:
    @pytest.mark.parametrize(
        ("name", "expected"),
        [
            ("a.csv | raw_data", SpectrumLabel("a.csv", "raw_data", ROLE_TITRATION)),
            ("Background | buffer", SpectrumLabel("buffer", "raw_data", ROLE_BUFFER)),
            (
                "Background | sol_A_buffer_subtracted",
                SpectrumLabel("sol_A", "buffer_subtracted", ROLE_SOL_A),
            ),
            (
                "Background | sol_B_buffer_subtracted_delta_epsilon",
                SpectrumLabel(
                    "sol_B", "buffer_subtracted_delta_epsilon", ROLE_SOL_B
                ),
            ),
            ("Background | other", SpectrumLabel("Background", "other", ROLE_TITRATION)),
            ("no_stage", SpectrumLabel("no_stage", "", ROLE_TITRATION)),
        ],
    )
    def test_names_without_label_metas_are_split(self, name, expected):
        table = spectra_table([name], np.zeros((1, 4)))
        assert spectrum_labels(table) == [expected]

    def test_label_name(self):
        assert SpectrumLabel("s", "t", ROLE_SOL_A).name == "s | t"

    def labelled(self):
        # Names that would be misread if they were split: a sample whose name
        # contains the separator, and a data sample called "Background".
        labels = [
            SpectrumLabel("sol_A", STAGE_RAW, ROLE_SOL_A),
            SpectrumLabel("odd | name.csv", STAGE_RAW, ROLE_TITRATION),
            SpectrumLabel("Background", "sol_A", ROLE_TITRATION),
        ]
        template = spectra_table(["x"], np.zeros((1, 4)))
        return build_spectra_table(
            template, np.zeros((3, 4)), labels, ["degree"] * 3
        ), labels

    def test_label_metas_win_over_the_name(self):
        table, labels = self.labelled()
        assert spectrum_labels(table) == labels
        assert [n for _, n in matching_spectra(table, role=ROLE_TITRATION)] == [
            "odd | name.csv | raw_data",
            "Background | sol_A",
        ]
        assert stages(table, ROLE_SOL_A) == [STAGE_RAW]
        assert stages(table, ROLE_TITRATION) == [STAGE_RAW, "sol_A"]

    def test_metas_are_written(self):
        table, _ = self.labelled()
        metas = [v.name for v in table.domain.metas]
        assert metas == ["Spectrum", "Sample", "Stage", "Role", "Unit"]


class TestUnits:
    def test_unit_strings(self):
        assert unit_string("nanometer") == "nanometer"
        assert unit_symbol("nanometer") == "nm"
        assert unit_symbol(None) is None
        assert unit_symbol("") is None

    def test_table_unit(self, table):
        assert table_unit(table, "spectrum_unit") == "millidegree"
        assert table_unit(table, "absent") is None

    def test_spectrum_units_from_attribute(self, table):
        assert spectrum_units(table) == ["millidegree"] * 4

    def test_spectrum_units_meta_wins(self):
        units = ["millidegree", "degree", "millidegree", "degree"]
        table = spectra_table(NAMES, np.zeros((4, 4)), units=units)
        assert spectrum_units(table) == units

    def test_spectrum_units_unknown(self):
        table = spectra_table(NAMES, np.zeros((4, 4)), spectrum_unit=None)
        assert spectrum_units(table) == [None] * 4

    def test_shared_unit(self):
        assert shared_unit(["a", "a"]) == "a"
        assert shared_unit(["a", "b"]) is None
        assert shared_unit([]) is None
        assert shared_unit(iter(["a"])) == "a"

    def test_delta_a_constant(self):
        assert MDEG_PER_DELTA_A.to("millidegree").magnitude == 32980.0


class TestQuantities:
    def test_round_trip(self, table):
        table.attributes["c"] = quantity_string(Q_(15, "micromolar"))
        quantity = table_quantity(table, "c")
        assert quantity.to("molar").magnitude == pytest.approx(15e-6)

    def test_round_trip_numpy_scalar(self, table):
        table.attributes["c"] = quantity_string(Q_(np.float64(2.5), "centimeter"))
        assert table_quantity(table, "c").magnitude == 2.5

    def test_missing_key(self, table):
        with pytest.raises(KeyError):
            table_quantity(table, "absent")

    @pytest.mark.parametrize("text", ["garbage ???", "0 micromolar", "-1 cm", "nan cm"])
    def test_invalid(self, table, text):
        table.attributes["c"] = text
        with pytest.raises(ValueError):
            table_quantity(table, "c")


class TestBuildSpectraTable:
    def test_layout_and_attributes(self, table):
        X = np.ones((2, 4))
        out = build_spectra_table(table, X, ["x | s", "y | s"], ["degree", "degree"])
        assert spectrum_names(out) == ["x | s", "y | s"]
        assert spectrum_units(out) == ["degree", "degree"]
        np.testing.assert_array_equal(wavelengths(out), wavelengths(table))
        assert out.attributes["wavelength_unit"] == "nanometer"
        assert out.attributes is not table.attributes

    def test_attribute_override(self, table):
        out = build_spectra_table(
            table, np.ones((1, 4)), ["x | s"], ["degree"], {"k": 1}
        )
        assert dict(out.attributes) == {"k": 1}


class TestPreprocessingDetection:
    """Quasar preprocessing leaves a compute_value on the attributes."""

    @pytest.fixture
    def preprocess(self):
        return pytest.importorskip("orangecontrib.spectroscopy.preprocess")

    def test_raw_table_is_not_preprocessed(self, table):
        assert not has_spectroscopy_preprocessing(table)

    def test_foreign_table(self):
        assert not has_spectroscopy_preprocessing(Table("iris"))

    def test_table_without_attributes(self):
        empty = Table.from_numpy(Domain([]), np.zeros((2, 0)))
        assert not has_spectroscopy_preprocessing(empty)

    def test_baseline_is_detected(self, table, preprocess):
        assert has_spectroscopy_preprocessing(preprocess.LinearBaseline()(table))

    def test_normalisation_is_detected(self, table, preprocess):
        normalised = preprocess.Normalize(method=preprocess.Normalize.MinMax)(table)
        assert has_spectroscopy_preprocessing(normalised)

    def test_detected_after_a_later_non_quasar_transform(self, table, preprocess):
        preprocessed = preprocess.LinearBaseline()(table)
        assert has_spectroscopy_preprocessing(preprocessed[:2])

    def test_cut_alone_is_not_detected(self, table, preprocess):
        # Cut only selects columns, so it leaves no trace (documented limit).
        assert not has_spectroscopy_preprocessing(
            preprocess.Cut(lowlim=235, highlim=255)(table)
        )

    def test_cut_then_baseline_is_detected(self, table, preprocess):
        cut = preprocess.Cut(lowlim=235, highlim=255)(table)
        assert has_spectroscopy_preprocessing(preprocess.LinearBaseline()(cut))
