import numpy as np
import pytest
from helpers import spectra_table
from Orange.widgets.tests import base as orange_tests

pytest.importorskip("cdpro")

from orangecontrib.orangeCD.widgets.core.owsecondarystructure import (
    OWSecondaryStructure,
    find_delta_epsilon_rows,
    wavelengths_nm,
)
from orangecontrib.orangeCD.widgets.titration.utils import (
    DELTA_EPSILON_UNIT,
    spectrum_names,
    spectrum_units,
    unit_string,
    wavelengths,
)

DE_UNIT = unit_string(DELTA_EPSILON_UNIT)
WL = np.arange(240.0, 189.0, -1.0)


def helix_like(scale=1.0):
    """A smooth helix-like delta epsilon spectrum, 240-190 nm in 1 nm steps."""
    return scale * (
        -8 * np.exp(-(((WL - 222) / 8) ** 2))
        - 10 * np.exp(-(((WL - 208) / 8) ** 2))
        + 20 * np.exp(-(((WL - 192) / 6) ** 2))
    )


def de_table(
    names=("a | plus_sol_A_delta_epsilon", "b | plus_sol_A_delta_epsilon"), **kw
):
    X = np.vstack([helix_like(1 + 0.1 * i) for i in range(len(names))])
    return spectra_table(names, X, units=[DE_UNIT] * len(names), wavelengths=WL, **kw)


class TestDetection:
    def test_rows_found_by_unit_and_name(self):
        table = spectra_table(
            ["x | plus_sol_A", "x | plus_sol_A_delta_epsilon", "y | other"],
            np.zeros((3, 4)),
            units=["millidegree", DE_UNIT, DE_UNIT],
        )
        found = find_delta_epsilon_rows(table)
        assert found == [(1, "x | plus_sol_A_delta_epsilon"), (2, "y | other")]

    def test_name_suffix_without_unit(self):
        table = spectra_table(["x | a_delta_epsilon", "x | a"], np.zeros((2, 4)))
        assert find_delta_epsilon_rows(table) == [(0, "x | a_delta_epsilon")]

    def test_wavelengths_converted_to_nm(self):
        table = spectra_table(["a"], np.zeros((1, 4)), wavelength_unit="micrometer")
        np.testing.assert_allclose(
            wavelengths_nm(table), np.array(wavelengths(table)) * 1000
        )


class TestWidget(orange_tests.WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWSecondaryStructure)

    def fit(self, table, method="SELCON3"):
        self.widget.method = method
        self.send_signal(self.widget.Inputs.data, table, widget=self.widget)
        self.wait_until_finished(self.widget, timeout=120000)
        return self.get_output(self.widget.Outputs.structure, widget=self.widget)

    def test_no_data(self):
        self.send_signal(self.widget.Inputs.data, None, widget=self.widget)
        assert (
            self.get_output(self.widget.Outputs.structure, widget=self.widget) is None
        )

    def test_missing_spectrum_meta(self):
        from Orange.data import Table

        self.send_signal(self.widget.Inputs.data, Table("iris"), widget=self.widget)
        assert self.widget.Error.missing_wavelength.is_shown()

    def test_no_delta_epsilon_rows(self):
        table = spectra_table(["a | raw_data"], np.zeros((1, 4)))
        self.send_signal(self.widget.Inputs.data, table, widget=self.widget)
        assert self.widget.Error.no_delta_epsilon.is_shown()

    def test_fit(self):
        for method in ("SELCON3", "CDSSTR", "CONTINLL"):
            with self.subTest(method=method):
                out = self.fit(de_table(), method)
                assert not self.widget.Error.active
                assert len(out) == 2
                assert spectrum_names(out) == ["a | plus_sol_A", "b | plus_sol_A"]
                assert not np.isnan(out.get_column("rmsd")).any()
                total = sum(
                    out.get_column(v)
                    for v in out.domain.attributes
                    if v.name.startswith("fraction_")
                )
                np.testing.assert_allclose(total, 1.0, atol=0.15)

    def test_fitted_spectra_layout(self):
        self.fit(de_table())
        spectra = self.get_output(self.widget.Outputs.spectra, widget=self.widget)
        assert spectrum_names(spectra) == [
            "a | plus_sol_A | measured",
            "a | plus_sol_A | calculated",
            "b | plus_sol_A | measured",
            "b | plus_sol_A | calculated",
        ]
        assert spectrum_units(spectra) == [DE_UNIT] * 4
        np.testing.assert_allclose(wavelengths(spectra), WL)
        assert spectra.attributes["wavelength_unit"] == "nanometer"
        np.testing.assert_allclose(spectra.X[0], helix_like(), atol=1e-6)

    def test_data_outside_reference_range_explains_why(self):
        wl = np.arange(300.0, 259.0, -1.0)
        table = spectra_table(
            ["a | s_delta_epsilon"],
            np.zeros((1, wl.size)),
            units=[DE_UNIT],
            wavelengths=wl,
        )
        out = self.fit(table)
        assert self.widget.Error.none_fitted.is_shown()
        assert "no overlap" in str(self.widget.Error.none_fitted)
        assert out is not None  # rows still reported in the output table

    def test_unfittable_spectrum_is_reported(self):
        names = ("a | s_delta_epsilon", "b | s_delta_epsilon")
        X = np.vstack([helix_like(), helix_like()])
        X[1, 5:7] = np.nan  # leaves a gap in the 1 nm grid
        table = spectra_table(names, X, units=[DE_UNIT] * 2, wavelengths=WL)
        out = self.fit(table)
        status = [
            str(m) for m in out.metas[:, out.domain.metas.index(out.domain["Status"])]
        ]
        assert status[0] == "ok"
        assert status[1] == "error"
        assert self.widget.Warning.some_failed.is_shown()
        assert "every nanometre" in str(self.widget.Warning.some_failed)
