#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Estimate protein secondary structure from delta epsilon CD spectra (CDPro)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.utils.concurrent import ConcurrentWidgetMixin, TaskState
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from cdpro import REFSETS
from cdpro import cdsstr, continll, selcon3
from cdpro.io.signal import to_dataset

from . import Q_

DELTA_EPSILON_SUFFIX = "delta_epsilon"
DELTA_EPSILON_UNIT = "liter / mole / centimeter"
AUTOMATIC = "Automatic (from wavelength range)"

METHODS = {
    "SELCON3": selcon3,
    "CDSSTR": cdsstr,
    "CONTINLL": continll,
}
REFERENCE_SETS = [AUTOMATIC] + [
    f"{info.name.strip()} - {info.description}" for info in REFSETS.values()
]
# Value passed to cdpro for each entry of REFERENCE_SETS.
REFERENCE_SET_VALUES = [None] + list(REFSETS)


# --------------------------------------------------------------------------
# Locating the wavelength and delta epsilon columns
# --------------------------------------------------------------------------

def find_wavelength_variable(data: Table) -> ContinuousVariable | None:
    """The continuous variable holding wavelengths (nm), if any."""
    candidates = [
        var for var in (*data.domain.metas, *data.domain.attributes)
        if isinstance(var, ContinuousVariable)
    ]
    for var in candidates:
        if var.name.strip().lower() in ("wavelength", "wavelength_nm", "wavelength (nm)"):
            return var
    for var in candidates:
        if var.name.strip().lower().startswith("wavelength"):
            return var
    return None


def _is_delta_epsilon(var: ContinuousVariable) -> bool:
    name = var.name.lower().replace(" ", "_")
    if name.endswith(DELTA_EPSILON_SUFFIX) or "delta_epsilon" in name or "Δε" in var.name:
        return True
    unit = var.attributes.get("unit")
    return bool(unit) and unit == str(Q_(1, DELTA_EPSILON_UNIT).units)


def find_delta_epsilon_variables(
    data: Table, wavelength: ContinuousVariable
) -> list[ContinuousVariable]:
    """Continuous variables that look like delta epsilon spectra."""
    return [
        var for var in (*data.domain.attributes, *data.domain.metas)
        if isinstance(var, ContinuousVariable)
        and var is not wavelength
        and _is_delta_epsilon(var)
    ]


def spectrum_label(variable: ContinuousVariable) -> str:
    name = variable.name
    if name.endswith(f"_{DELTA_EPSILON_SUFFIX}"):
        name = name[: -len(DELTA_EPSILON_SUFFIX) - 1]
    return name


# --------------------------------------------------------------------------
# Fitting
# --------------------------------------------------------------------------

@dataclass
class FitSettings:
    method: str
    basis: int | None
    wl_min: float | None
    wl_max: float | None
    seed: int


@dataclass
class SpectrumResult:
    label: str
    status: str = "error"
    message: str = ""
    reference_set: str = ""
    fractions: dict[str, float] = field(default_factory=dict)
    totals: dict[str, float] = field(default_factory=dict)
    rmsd: float = np.nan
    nrmsd: float = np.nan
    wavelength: np.ndarray = field(default_factory=lambda: np.empty(0))
    measured: np.ndarray = field(default_factory=lambda: np.empty(0))
    calculated: np.ndarray = field(default_factory=lambda: np.empty(0))


def fit_spectrum(
    label: str, wavelength: np.ndarray, delta_epsilon: np.ndarray, settings: FitSettings
) -> SpectrumResult:
    """Run the chosen cdpro method on one spectrum."""
    result = SpectrumResult(label)
    keep = np.isfinite(wavelength) & np.isfinite(delta_epsilon)
    try:
        dataset = to_dataset(wavelength[keep], delta_epsilon[keep], title=label)
        kwargs = {"basis": settings.basis, "wl_min": settings.wl_min, "wl_max": settings.wl_max}
        if settings.method == "CDSSTR":
            kwargs["seed"] = settings.seed
        doc = METHODS[settings.method].to_json(dataset, **kwargs)
        entry = doc["results"][0]
    except Exception as exc:  # cdpro raises ValueError for unusable spectra
        result.message = str(exc)
        return result

    result.status = entry["status"]
    result.message = entry.get("message") or ""
    result.reference_set = entry["reference_set"]["name"]
    if entry["status"] != "ok":
        return result
    result.fractions = dict(entry["secondary_structure"])
    result.totals = dict(entry.get("totals") or {})
    result.rmsd = entry["fit"]["rmsd"]
    result.nrmsd = entry["fit"]["nrmsd"]
    spectra = entry["spectra"]
    result.wavelength = np.asarray(spectra["wavelength_nm"], dtype=float)
    result.measured = np.asarray(spectra["measured"], dtype=float)
    result.calculated = np.asarray(spectra["calculated"], dtype=float)
    return result


def fit_all(
    wavelength: np.ndarray,
    spectra: list[tuple[str, np.ndarray]],
    settings: FitSettings,
    state: TaskState,
) -> list[SpectrumResult]:
    results = []
    for i, (label, values) in enumerate(spectra):
        if state.is_interruption_requested():
            break
        state.set_status(f"Fitting {label}")
        results.append(fit_spectrum(label, wavelength, values, settings))
        state.set_progress_value(100 * (i + 1) / len(spectra))
    return results


# --------------------------------------------------------------------------
# Building output tables
# --------------------------------------------------------------------------

def _ordered_union(groups: list[dict]) -> list[str]:
    keys: list[str] = []
    for group in groups:
        for key in group:
            if key not in keys:
                keys.append(key)
    return keys


def structure_table(results: list[SpectrumResult], name: str = "") -> Table:
    fraction_keys = _ordered_union([r.fractions for r in results])
    total_keys = _ordered_union([r.totals for r in results])
    attributes = (
        [ContinuousVariable(f"fraction_{key}") for key in fraction_keys]
        + [ContinuousVariable(f"total_{key}") for key in total_keys]
        + [ContinuousVariable("rmsd"), ContinuousVariable("nrmsd")]
    )
    metas = [
        StringVariable("Spectrum"), StringVariable("Reference set"),
        StringVariable("Status"), StringVariable("Message"),
    ]
    X = np.array([
        [r.fractions.get(k, np.nan) for k in fraction_keys]
        + [r.totals.get(k, np.nan) for k in total_keys]
        + [r.rmsd, r.nrmsd]
        for r in results
    ], dtype=float).reshape(len(results), len(attributes))
    M = np.array(
        [[r.label, r.reference_set, r.status, r.message] for r in results], dtype=object
    ).reshape(len(results), len(metas))
    table = Table.from_numpy(Domain(attributes, metas=metas), X, metas=M)
    table.name = name or "Secondary structure"
    return table


def fitted_spectra_table(results: list[SpectrumResult], name: str = "") -> Table | None:
    """Measured and calculated spectra of the fits on a common wavelength grid."""
    fitted = [r for r in results if r.wavelength.size]
    if not fitted:
        return None
    grid = np.unique(np.concatenate([r.wavelength for r in fitted]))[::-1]
    columns, variables = [], []
    for r in fitted:
        index = np.searchsorted(-grid, -r.wavelength)
        for kind, values in (("measured", r.measured), ("calculated", r.calculated)):
            column = np.full(grid.size, np.nan)
            column[index] = values
            columns.append(column)
            variable = ContinuousVariable(f"{r.label} | {kind}")
            variable.attributes["unit"] = str(Q_(1, DELTA_EPSILON_UNIT).units)
            variables.append(variable)
    wavelength = ContinuousVariable("Wavelength")
    table = Table.from_numpy(
        Domain(variables, metas=[wavelength]),
        np.column_stack(columns),
        metas=grid.reshape(-1, 1),
    )
    table.name = name or "Fitted spectra"
    return table


# --------------------------------------------------------------------------
# Widget
# --------------------------------------------------------------------------

class OWSecondaryStructure(OWWidget, ConcurrentWidgetMixin):
    name = "Secondary Structure"
    description = (
        "Estimate protein secondary structure fractions from delta epsilon "
        "CD spectra using CDPro methods."
    )
    icon = "icons/SecondaryStructure.svg"
    priority = 40
    want_main_area = False
    resizing_enabled = False

    class Inputs:
        data = Input("CD and Delta Epsilon Spectra", Table)

    class Outputs:
        structure = Output("Secondary Structure", Table)
        spectra = Output("Fitted Spectra", Table)

    method = Setting("CDSSTR")
    reference_set_index = Setting(0)
    use_wl_min = Setting(False)
    wl_min = Setting(190.0)
    use_wl_max = Setting(False)
    wl_max = Setting(240.0)
    seed = Setting(0)
    auto_commit = Setting(False)

    class Error(OWWidget.Error):
        missing_wavelength = Msg("Input does not contain a wavelength column.")
        no_delta_epsilon = Msg(
            "No delta epsilon columns found (names ending in '_delta_epsilon' "
            "or carrying the delta epsilon unit)."
        )
        fit_failed = Msg("{}")

    class Warning(OWWidget.Warning):
        some_failed = Msg("{} of {} spectra could not be fitted; see the Status column.")

    def __init__(self) -> None:
        OWWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)
        if self.method not in METHODS:
            self.method = "CDSSTR"
        if not 0 <= self.reference_set_index < len(REFERENCE_SETS):
            self.reference_set_index = 0

        self.data: Table | None = None
        self._wavelength: np.ndarray | None = None
        self._spectra: list[tuple[str, np.ndarray]] = []
        self._build_controls()

    def _build_controls(self) -> None:
        box = gui.widgetBox(self.controlArea, "Method")
        gui.comboBox(
            box, self, "method", label="Fitting method", items=list(METHODS),
            sendSelectedValue=True, orientation="horizontal",
            callback=self._settings_changed,
        )
        gui.comboBox(
            box, self, "reference_set_index", label="Reference set",
            items=REFERENCE_SETS, orientation="horizontal",
            callback=self._settings_changed,
        )
        self.seed_spin = gui.spin(
            box, self, "seed", 0, 2**31 - 1, label="CDSSTR random seed",
            orientation="horizontal", callback=self._settings_changed,
        )

        limits = gui.widgetBox(self.controlArea, "Wavelength limits (nm)")
        gui.checkBox(limits, self, "use_wl_min", "Shortest wavelength",
                     callback=self._settings_changed)
        gui.doubleSpin(limits, self, "wl_min", 100, 400, step=1, decimals=1,
                       callback=self._settings_changed)
        gui.checkBox(limits, self, "use_wl_max", "Longest wavelength",
                     callback=self._settings_changed)
        gui.doubleSpin(limits, self, "wl_max", 100, 400, step=1, decimals=1,
                       callback=self._settings_changed)

        self.info_label = gui.widgetLabel(self.controlArea, "No data.")
        self.info_label.setWordWrap(True)
        note = gui.widgetLabel(
            self.controlArea,
            "Each delta epsilon column is fitted as a separate spectrum. "
            "Data must cover every whole nanometre in its range.",
        )
        note.setWordWrap(True)
        gui.auto_commit(
            self.buttonsArea, self, "auto_commit", "Fit", commit=self.commit
        )
        self._update_seed_enabled()

    def _update_seed_enabled(self) -> None:
        self.seed_spin.setEnabled(self.method == "CDSSTR")

    def _settings_changed(self) -> None:
        self._update_seed_enabled()
        self.commit.deferred()

    @Inputs.data
    def set_data(self, data: Table | None) -> None:
        self.data = data
        self._wavelength = None
        self._spectra = []
        self.Error.clear()
        self.Warning.clear()
        self.cancel()
        if data is not None:
            self._locate_columns(data)
        self.info_label.setText(
            f"{len(self._spectra)} delta epsilon spectra found."
            if self._spectra else "No spectra to fit."
        )
        self.commit.now()

    def _locate_columns(self, data: Table) -> None:
        wavelength = find_wavelength_variable(data)
        if wavelength is None:
            self.Error.missing_wavelength()
            return
        variables = find_delta_epsilon_variables(data, wavelength)
        if not variables:
            self.Error.no_delta_epsilon()
            return
        self._wavelength = np.asarray(data.get_column(wavelength), dtype=float)
        self._spectra = [
            (spectrum_label(var), np.asarray(data.get_column(var), dtype=float))
            for var in variables
        ]

    def _settings(self) -> FitSettings:
        return FitSettings(
            method=self.method,
            basis=REFERENCE_SET_VALUES[self.reference_set_index],
            wl_min=self.wl_min if self.use_wl_min else None,
            wl_max=self.wl_max if self.use_wl_max else None,
            seed=int(self.seed),
        )

    @gui.deferred
    def commit(self) -> None:
        self.Error.fit_failed.clear()
        self.Warning.some_failed.clear()
        self.cancel()
        if self._wavelength is None or not self._spectra:
            self._send(None, None)
            return
        self.start(fit_all, self._wavelength, self._spectra, self._settings())

    def on_partial_result(self, result) -> None:
        pass

    def on_done(self, results: list[SpectrumResult]) -> None:
        name = self.data.name if self.data is not None and self.data.name else ""
        failed = sum(r.status != "ok" for r in results)
        if failed:
            self.Warning.some_failed(failed, len(results))
        self._send(
            structure_table(results, f"{name} - secondary structure" if name else ""),
            fitted_spectra_table(results, f"{name} - fitted spectra" if name else ""),
        )

    def on_exception(self, ex: Exception) -> None:
        self.Error.fit_failed(str(ex))
        self._send(None, None)

    def _send(self, structure: Table | None, spectra: Table | None) -> None:
        self.Outputs.structure.send(structure)
        self.Outputs.spectra.send(spectra)

    def onDeleteWidget(self) -> None:
        self.shutdown()
        super().onDeleteWidget()


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWSecondaryStructure).run()
