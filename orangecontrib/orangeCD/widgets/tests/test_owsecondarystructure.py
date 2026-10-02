import numpy as np
import pytest

pytest.importorskip("cdpro")

from cdpro.io.input_file import parse_input
from Orange.data import ContinuousVariable, Domain, Table
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.orangeCD.widgets.titration.owsecondarystructure import (
    OWSecondaryStructure,
    find_delta_epsilon_variables,
    find_wavelength_variable,
)


def _spectrum():
    """A smooth helix-like delta epsilon spectrum, 240-190 nm in 1 nm steps."""
    wl = np.arange(240.0, 189.0, -1.0)
    de = (
        -8 * np.exp(-((wl - 222) / 8) ** 2)
        - 10 * np.exp(-((wl - 208) / 8) ** 2)
        + 20 * np.exp(-((wl - 192) / 6) ** 2)
    )
    return wl, de


def _table(columns=("a | plus_sol_A_delta_epsilon", "b | plus_sol_A_delta_epsilon")):
    wl, de = _spectrum()
    variables = [ContinuousVariable(name) for name in columns]
    X = np.column_stack([de * (1 + 0.1 * i) for i in range(len(columns))])
    return Table.from_numpy(
        Domain(variables, metas=[ContinuousVariable("Wavelength")]), X,
        metas=wl.reshape(-1, 1),
    )


def test_column_detection():
    table = _table(("x | plus_sol_A", "x | plus_sol_A_delta_epsilon"))
    wavelength = find_wavelength_variable(table)
    assert wavelength.name == "Wavelength"
    found = find_delta_epsilon_variables(table, wavelength)
    assert [v.name for v in found] == ["x | plus_sol_A_delta_epsilon"]


class TestOWSecondaryStructure(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWSecondaryStructure)

    def test_no_data(self):
        self.send_signal(self.widget.Inputs.data, None)
        assert self.get_output(self.widget.Outputs.structure) is None

    def test_missing_columns(self):
        wl, de = _spectrum()
        table = Table.from_numpy(Domain([ContinuousVariable("other")]), de.reshape(-1, 1))
        self.send_signal(self.widget.Inputs.data, table)
        assert self.widget.Error.missing_wavelength.is_shown()

    @pytest.mark.parametrize("method", ["SELCON3", "CDSSTR", "CONTINLL"])
    def test_fit(self, method):
        self.widget.method = method
        self.send_signal(self.widget.Inputs.data, _table())
        self.wait_until_finished(self.widget, timeout=120000)
        out = self.get_output(self.widget.Outputs.structure)
        assert len(out) == 2
        assert [str(m[0]) for m in out.metas] == ["a | plus_sol_A", "b | plus_sol_A"]
        assert out.domain["rmsd"] is not None
        spectra = self.get_output(self.widget.Outputs.spectra)
        assert "a | plus_sol_A | calculated" in [v.name for v in spectra.domain.attributes]
