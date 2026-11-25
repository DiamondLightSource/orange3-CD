#!/usr/bin/env python3
"""Apply the titration reference corrections to loaded CD spectra."""

from __future__ import annotations

import numpy as np
from AnyQt.QtCore import Qt
from Orange.data import Table
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from ...units import Q_
from .utils import (
    ROLE_BUFFER,
    ROLE_SOL_A,
    ROLE_SOL_B,
    ROLE_TITRATION,
    STAGE_BUFFER_SUBTRACTED,
    STAGE_RAW,
    STAGE_SOL_A_ADDED,
    STAGE_SOL_A_SUBTRACTED,
    STAGE_SOL_B_SUBTRACTED,
    WAVELENGTH_UNIT_KEY,
    InvalidWavelength,
    SpectraError,
    SpectrumLabel,
    build_spectra_table,
    has_spectroscopy_preprocessing,
    spectrum_labels,
    spectrum_points,
    spectrum_units,
    table_unit,
    titration_rows,
    unit_symbol,
)
from .utils import (
    wavelengths as spectra_wavelengths,
)

DATA_STAGE = STAGE_RAW
DILUTION_COLUMN = "dilution_factor"
RATIO_COLUMN = "normalised_molar_ratio"

# Stages written for every titration data spectrum, in order.
DATA_STAGES = (
    STAGE_BUFFER_SUBTRACTED,
    STAGE_SOL_A_SUBTRACTED,
    STAGE_SOL_B_SUBTRACTED,
    STAGE_SOL_A_ADDED,
)


def correct_spectra(
    raw: np.ndarray,
    buffer: np.ndarray,
    sol_a: np.ndarray,
    sol_b: np.ndarray | None,
    dilution: np.ndarray,
    ratio: np.ndarray,
) -> dict[str, np.ndarray]:
    """Apply the titration corrections; every spectrum shares one unit.

    ``raw`` is ``(n_points, n_wavelengths)``; the references are
    ``(n_wavelengths,)`` and ``dilution`` / ``ratio`` are ``(n_points,)``.
    Returns the new reference rows (``sol_A_buffer_subtracted``, and
    ``sol_B_buffer_subtracted`` when Solution B is given) and one
    ``(n_points, n_wavelengths)`` array per entry of ``DATA_STAGES``.

    1. buffer subtraction, scaled by the dilution factor of the point;
    2. subtract buffer-subtracted Solution A;
    3. subtract buffer-subtracted Solution B scaled by the normalised molar
       ratio (skipped without a Solution B spectrum);
    4. add the buffer-subtracted Solution A back.

    The zero-level offset is not applied here; do it in spectral
    preprocessing beforehand.
    """
    sol_a_bs = sol_a - buffer
    result = {"sol_A_buffer_subtracted": sol_a_bs}
    buffer_subtraction = (raw - buffer) * dilution[:, None]
    sol_a_subtraction = buffer_subtraction - sol_a_bs
    subtract_frac_sol_b = sol_a_subtraction
    if sol_b is not None:
        sol_b_bs = sol_b - buffer
        result["sol_B_buffer_subtracted"] = sol_b_bs
        subtract_frac_sol_b = sol_a_subtraction - sol_b_bs * ratio[:, None]
    result.update(
        {
            STAGE_BUFFER_SUBTRACTED: buffer_subtraction,
            STAGE_SOL_A_SUBTRACTED: sol_a_subtraction,
            STAGE_SOL_B_SUBTRACTED: subtract_frac_sol_b,
            STAGE_SOL_A_ADDED: subtract_frac_sol_b + sol_a_bs,
        }
    )
    return result


DEFAULT_LAZY_WINDOW_PERCENT = 10
MIN_WINDOW_POINTS = 3

# Zero-level baseline subtracted before the corrections.
BASELINE_NONE, BASELINE_AUTO, BASELINE_RANGE = range(3)
# Fixed-range defaults, in the wavelength unit of the input (normally nm):
# a stretch above the far-UV protein bands where CD should be zero.
DEFAULT_RANGE_START = 250.0
DEFAULT_RANGE_END = 260.0


def find_flat_region(X: np.ndarray, window: int) -> tuple[int, int]:
    """Find the flattest stretch of the wavelength axis.

    Every window of ``window`` consecutive points is scored by the spread of
    each spectrum inside it, relative to that spectrum's spread over the
    whole axis, so that each spectrum counts equally whatever its size.
    Returns the half-open ``(start, stop)`` index range of the best window.
    """
    X = np.asarray(X, dtype=float)
    n_points = X.shape[1]
    window = max(1, min(window, n_points))
    overall = np.nanstd(X, axis=1)
    usable = overall > 0
    if not usable.any():
        return 0, window
    windows = np.lib.stride_tricks.sliding_window_view(X[usable], window, axis=1)
    score = np.nanmean(np.nanstd(windows, axis=2) / overall[usable, None], axis=0)
    start = int(np.nanargmin(score))
    return start, start + window


def auto_baseline(
    X: np.ndarray, window_percent: float = DEFAULT_LAZY_WINDOW_PERCENT
) -> tuple[np.ndarray, np.ndarray, tuple[int, int]]:
    """Subtract each spectrum's mean over the automatically found flat region.

    Returns ``(corrected X, per-spectrum offsets, (start, stop))``.
    """
    X = np.asarray(X, dtype=float)
    window = max(MIN_WINDOW_POINTS, round(X.shape[1] * window_percent / 100))
    start, stop = find_flat_region(X, window)
    offsets = np.nanmean(X[:, start:stop], axis=1)
    return X - offsets[:, None], offsets, (start, stop)


def range_baseline(
    X: np.ndarray, wavelengths: np.ndarray, start: float, end: float
) -> tuple[np.ndarray, np.ndarray, tuple[float, float]]:
    """Subtract each spectrum's mean between two wavelengths (inclusive).

    ``start`` and ``end`` may be given in either order. Returns
    ``(corrected X, per-spectrum offsets, (lowest, highest) wavelength used)``
    and raises ``ValueError`` when no wavelength falls inside the range.
    """
    X = np.asarray(X, dtype=float)
    wavelengths = np.asarray(wavelengths, dtype=float)
    low, high = sorted((start, end))
    inside = (wavelengths >= low) & (wavelengths <= high)
    if not inside.any():
        raise ValueError(
            f"No wavelengths between {low:g} and {high:g}; the spectra cover "
            f"{wavelengths.min():g}\u2013{wavelengths.max():g}"
        )
    offsets = np.nanmean(X[:, inside], axis=1)
    used = wavelengths[inside]
    return X - offsets[:, None], offsets, (float(used.min()), float(used.max()))


class OWCDDataCorrection(OWWidget):
    name = "CD Data Correction"
    description = (
        "Apply buffer, Solution A and Solution B corrections to preprocessed "
        "CD titration spectra."
    )
    icon = "icons/CDDataCorrection.svg"
    priority = 25
    want_main_area = False
    resizing_enabled = False

    settings_version = 2
    baseline_mode = Setting(BASELINE_NONE)
    auto_baseline_window_percent = Setting(DEFAULT_LAZY_WINDOW_PERCENT)
    min_wavelength = Setting(DEFAULT_RANGE_START)
    max_wavelength = Setting(DEFAULT_RANGE_END)

    class Inputs:
        data = Input("CD Data", Table)
        titration = Input("Titration Table", Table)

    class Outputs:
        data = Output("Corrected CD Data", Table)

    class Error(OWWidget.Error):
        missing_spectrum_meta = Msg(
            "Input does not contain a 'Spectrum' string meta naming the rows."
        )
        invalid_wavelength = Msg("Input attribute names are not all valid wavelengths.")
        missing_background = Msg("No raw '{}' reference spectrum was found.")
        no_data_spectra = Msg(f"No '{DATA_STAGE}' spectra were found.")
        invalid_titration = Msg("{}")
        incompatible_units = Msg("{}")
        empty_baseline_range = Msg("{}")

    class Warning(OWWidget.Warning):
        no_titration = Msg("Connect a Titration Table to correct the spectra.")
        no_dilution_factor = Msg(
            f"The Titration Table has no '{DILUTION_COLUMN}' column: a dilution "
            "factor of 1 was assumed for every point (correct only for a "
            "fixed-volume titration)."
        )
        no_preprocessing = Msg(
            "The input does not appear to have been through Preprocess Spectra,\n"
            "so its baseline has not been subtracted. Use Preprocess Spectra,\n"
            "or choose a zero-level baseline here."
        )
        double_baseline = Msg(
            "The input has already been preprocessed, and the zero-level "
            "baseline chosen here subtracts a further one."
        )
        no_sol_b = Msg(
            "No raw 'sol_B' reference spectrum: Solution B subtraction skipped."
        )

    def __init__(self) -> None:
        super().__init__()
        self.data: Table | None = None
        self.titration: Table | None = None
        box = gui.widgetBox(self.controlArea, "Corrections")
        note = gui.widgetLabel(
            box,
            "Applied to every titration 'raw_data' spectrum using the raw "
            "'buffer', 'sol_A' and 'sol_B' reference spectra and the dilution factor and "
            "normalised molar ratio of the Titration Table:\n"
            "1. buffer subtraction x dilution factor (1 if the Titration Table has none)\n"
            "2. subtract Solution A\n"
            "3. subtract Solution B x molar ratio\n"
            "4. add Solution A back\n"
            "The zero-level offset is not applied; use spectral preprocessing "
            "first, or choose a zero-level baseline below.",
        )
        note.setWordWrap(True)

        baseline_box = gui.radioButtons(
            self.controlArea,
            self,
            "baseline_mode",
            box="Zero-level baseline",
            callback=self._baseline_mode_changed,
        )
        gui.appendRadioButton(baseline_box, "None (input already preprocessed)")
        gui.appendRadioButton(
            baseline_box,
            "Flattest region",
            tooltip="Find the flattest region of the spectra and subtract each "
            "spectrum's mean over it before the corrections.",
        )
        self.auto_baseline_gui = gui.indentedBox(baseline_box)
        gui.spin(
            self.auto_baseline_gui,
            self,
            "auto_baseline_window_percent",
            2,
            50,
            label="Width (% of wavelength range)",
            orientation=Qt.Horizontal,
            callback=self.commit,
        )
        gui.appendRadioButton(
            baseline_box,
            "Fixed wavelength range",
            tooltip="Subtract each spectrum's mean between the two wavelengths "
            "before the corrections.",
        )
        self.range_baseline_gui = gui.indentedBox(baseline_box)
        self.range_spins = [
            gui.doubleSpin(
                self.range_baseline_gui,
                self,
                setting,
                0,
                10_000,
                0.5,
                label=label,
                orientation=Qt.Horizontal,
                callback=self.commit,
            )
            for setting, label in (
                ("min_wavelength", "From"),
                ("max_wavelength", "To"),
            )
        ]
        self.baseline_label = gui.widgetLabel(baseline_box, "")
        self.baseline_label.setWordWrap(True)
        self._update_baseline_controls()

    @classmethod
    def migrate_settings(cls, settings, version):
        if version is None or version < 2:
            # Version 1 had independent "lazy process" / "simple subtraction"
            # check boxes; they are now one choice.
            lazy = settings.pop("lazy_process", False)
            auto = settings.pop("auto_baseline_setting", False)
            simple = settings.pop("simple_subtraction_setting", False)
            if simple:
                settings["baseline_mode"] = BASELINE_RANGE
            elif lazy or auto:
                settings["baseline_mode"] = BASELINE_AUTO
            if "lazy_window_percent" in settings:
                settings.setdefault(
                    "auto_baseline_window_percent",
                    settings.pop("lazy_window_percent"),
                )
            # Placeholder range of the first fixed-range version.
            if (
                settings.get("min_wavelength") == 10
                and settings.get("max_wavelength") == 20
            ):
                del settings["min_wavelength"], settings["max_wavelength"]

    @Inputs.data
    def set_data(self, data: Table | None) -> None:
        self.data = data
        unit = None if data is None else table_unit(data, WAVELENGTH_UNIT_KEY)
        suffix = f" {unit_symbol(unit)}" if unit else ""
        for spin in self.range_spins:
            spin.setSuffix(suffix)
        self.commit()

    @Inputs.titration
    def set_titration(self, table: Table | None) -> None:
        self.titration = table
        self.commit()

    def _titration_column(self, name: str) -> np.ndarray:
        try:
            return np.asarray(self.titration.get_column(name), dtype=float)
        except (KeyError, ValueError) as exc:
            raise ValueError(
                f"Titration Table is missing a numeric '{name}' column"
            ) from exc

    def _update_baseline_controls(self) -> None:
        self.auto_baseline_gui.setEnabled(self.baseline_mode == BASELINE_AUTO)
        self.range_baseline_gui.setEnabled(self.baseline_mode == BASELINE_RANGE)

    def _baseline_mode_changed(self) -> None:
        self._update_baseline_controls()
        self.commit()

    def commit(self) -> None:
        self.Error.clear()
        self.Warning.clear()
        self.Outputs.data.send(self._correct())

    def _correct(self) -> Table | None:
        self.baseline_label.setText("")
        if self.data is None:
            return None
        preprocessed = has_spectroscopy_preprocessing(self.data)
        baseline = self.baseline_mode != BASELINE_NONE
        if baseline and preprocessed:
            self.Warning.double_baseline()
        elif not baseline and not preprocessed:
            self.Warning.no_preprocessing()
        if self.titration is None:
            self.Warning.no_titration()
            return None
        try:
            spectra_wavelengths(self.data)
            labels = spectrum_labels(self.data)
        except InvalidWavelength:
            self.Error.invalid_wavelength()
            return None
        except SpectraError:
            self.Error.missing_spectrum_meta()
            return None

        def reference_row(role: str) -> int | None:
            """The raw spectrum of a reference, if there is one."""
            for index, label in enumerate(labels):
                if label.role == role and label.stage == STAGE_RAW:
                    return index
            return None

        buffer_row = reference_row(ROLE_BUFFER)
        sol_a_row = reference_row(ROLE_SOL_A)
        sol_b_row = reference_row(ROLE_SOL_B)
        for role, row in ((ROLE_BUFFER, buffer_row), (ROLE_SOL_A, sol_a_row)):
            if row is None:
                self.Error.missing_background(role)
                return None
        if sol_b_row is None:
            self.Warning.no_sol_b()

        data_rows = [
            index
            for index, label in enumerate(labels)
            if label.role == ROLE_TITRATION and label.stage == DATA_STAGE
        ]
        if not data_rows:
            self.Error.no_data_spectra()
            return None

        points = spectrum_points(self.data)
        try:
            ratio = self._titration_column(RATIO_COLUMN)
            if DILUTION_COLUMN in self.titration.domain:
                dilution = self._titration_column(DILUTION_COLUMN)
            else:
                # Fixed-volume titrations do not dilute the cell.
                self.Warning.no_dilution_factor()
                dilution = np.ones(len(ratio))
            if points is None:
                # Unnumbered spectra are paired with the table rows in order.
                if len(ratio) != len(data_rows):
                    raise ValueError(
                        f"The Titration Table has {len(ratio)} rows but there "
                        f"are {len(data_rows)} '{DATA_STAGE}' spectra"
                    )
            else:
                rows = titration_rows(points[data_rows], len(ratio), DATA_STAGE)
                ratio, dilution = ratio[rows], dilution[rows]
        except ValueError as exc:
            self.Error.invalid_titration(str(exc))
            return None

        # Bring every spectrum to the unit of the first data spectrum.
        row_units = spectrum_units(self.data)
        unit = row_units[data_rows[0]]
        needed = [*data_rows, buffer_row, sol_a_row]
        if sol_b_row is not None:
            needed.append(sol_b_row)
        if not all(row_units[row] for row in needed):
            self.Error.incompatible_units(
                "The input has no spectrum unit "
                "(attribute 'spectrum_unit' or 'Unit' meta)."
            )
            return None

        X = self.data.X
        region = None
        if self.baseline_mode == BASELINE_AUTO:
            X, _, (start, stop) = auto_baseline(X, self.auto_baseline_window_percent)
            axis = spectra_wavelengths(self.data)[start:stop]
            region = (float(axis.min()), float(axis.max()))
            description = "the flat region"
        elif self.baseline_mode == BASELINE_RANGE:
            try:
                X, _, region = range_baseline(
                    X,
                    spectra_wavelengths(self.data),
                    self.min_wavelength,
                    self.max_wavelength,
                )
            except ValueError as exc:
                self.Error.empty_baseline_range(str(exc))
                return None
            description = "the range"
        if region is not None:
            symbol = unit_symbol(table_unit(self.data, WAVELENGTH_UNIT_KEY)) or ""
            self.baseline_label.setText(
                f"Baseline taken from {description} {region[0]:g}\u2013"
                f"{region[1]:g} {symbol}".rstrip()
            )

        def values(row: int) -> np.ndarray:
            return Q_(X[row], row_units[row]).to(unit).magnitude

        try:
            corrected = correct_spectra(
                np.vstack([values(row) for row in data_rows]),
                values(buffer_row),
                values(sol_a_row),
                None if sol_b_row is None else values(sol_b_row),
                dilution,
                ratio,
            )
        except TypeError as exc:  # pint dimensionality mismatch
            self.Error.incompatible_units(str(exc))
            return None

        samples = [labels[row].sample for row in data_rows]
        new_names, new_rows, new_points = [], [], []
        for role, key in (
            (ROLE_SOL_A, "sol_A_buffer_subtracted"),
            (ROLE_SOL_B, "sol_B_buffer_subtracted"),
        ):
            if key in corrected:
                new_names.append(SpectrumLabel(role, STAGE_BUFFER_SUBTRACTED, role))
                new_rows.append(corrected[key][None, :])
                new_points.append(np.nan)  # references belong to no point
        for stage in DATA_STAGES:
            new_names.extend(
                SpectrumLabel(sample, stage, ROLE_TITRATION) for sample in samples
            )
            new_rows.append(corrected[stage])
            if points is not None:
                new_points.extend(points[data_rows])

        # Keep the input rows so the original stages remain selectable.
        output = build_spectra_table(
            self.data,
            np.vstack((self.data.X, *new_rows)),
            [*labels, *new_names],
            [*row_units, *[unit] * len(new_names)],
            # Each corrected data spectrum keeps the point it came from.
            points=None if points is None else [*points, *new_points],
        )
        if region is not None:
            output.attributes["baseline_region"] = list(region)
        output.name = (
            f"{self.data.name} - corrected" if self.data.name else "Corrected CD data"
        )
        return output


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWCDDataCorrection).run()
