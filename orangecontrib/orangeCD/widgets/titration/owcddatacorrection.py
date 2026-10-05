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
    BACKGROUND_PREFIX,
    WAVELENGTH_UNIT_KEY,
    InvalidWavelength,
    SpectraError,
    build_spectra_table,
    has_spectroscopy_preprocessing,
    spectrum_names,
    spectrum_units,
    split_series_name,
    table_unit,
    unit_symbol,
)
from .utils import (
    wavelengths as spectra_wavelengths,
)

DATA_STAGE = "raw_data"
DILUTION_COLUMN = "dilution_factor"
RATIO_COLUMN = "normalised_molar_ratio"

# Stages written for every titration data spectrum, in order.
DATA_STAGES = (
    "buffer_subtraction",
    "sol_A_subtraction",
    "subtract_frac_sol_B",
    "plus_sol_A",
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
            "buffer_subtraction": buffer_subtraction,
            "sol_A_subtraction": sol_a_subtraction,
            "subtract_frac_sol_B": subtract_frac_sol_b,
            "plus_sol_A": subtract_frac_sol_b + sol_a_bs,
        }
    )
    return result


DEFAULT_LAZY_WINDOW_PERCENT = 10
MIN_WINDOW_POINTS = 3


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

    lazy_process = Setting(False)
    lazy_window_percent = Setting(DEFAULT_LAZY_WINDOW_PERCENT)

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
        missing_background = Msg("Background spectrum '{}' was not found.")
        no_data_spectra = Msg(f"No '{DATA_STAGE}' spectra were found.")
        invalid_titration = Msg("{}")
        incompatible_units = Msg("{}")

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
            "or press 'Lazy process' to subtract a baseline here."
        )
        double_baseline = Msg(
            "The input has already been preprocessed, and 'Lazy process' "
            "subtracts a further baseline."
        )
        no_sol_b = Msg(
            "No 'Background | sol_B' spectrum: Solution B subtraction skipped."
        )

    def __init__(self) -> None:
        super().__init__()
        self.data: Table | None = None
        self.titration: Table | None = None
        box = gui.widgetBox(self.controlArea, "Corrections")
        note = gui.widgetLabel(
            box,
            "Applied to every 'raw_data' spectrum using the 'Background | "
            "buffer', 'sol_A' and 'sol_B' spectra and the dilution factor and "
            "normalised molar ratio of the Titration Table:\n"
            "1. buffer subtraction x dilution factor (1 if the Titration Table has none)\n"
            "2. subtract Solution A\n"
            "3. subtract Solution B x molar ratio\n"
            "4. add Solution A back\n"
            "The zero-level offset is not applied; use spectral preprocessing "
            "first, or 'Lazy process' below.",
        )
        note.setWordWrap(True)

        lazy_box = gui.widgetBox(self.controlArea, "Lazy process")
        gui.checkBox(
            lazy_box,
            self,
            "lazy_process",
            "Lazy process",
            callback=self.commit,
            tooltip="Find the flattest region of the spectra and subtract each "
            "spectrum's mean over it before the corrections.",
        )
        gui.spin(
            lazy_box,
            self,
            "lazy_window_percent",
            2,
            50,
            label="Flat region width (% of wavelength range)",
            orientation=Qt.Horizontal,
            callback=self.commit,
        )
        self.lazy_label = gui.widgetLabel(lazy_box, "")
        self.lazy_label.setWordWrap(True)

    @Inputs.data
    def set_data(self, data: Table | None) -> None:
        self.data = data
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

    def commit(self) -> None:
        self.Error.clear()
        self.Warning.clear()
        self.Outputs.data.send(self._correct())

    def _correct(self) -> Table | None:
        self.lazy_label.setText("")
        if self.data is None:
            return None
        preprocessed = has_spectroscopy_preprocessing(self.data)
        if self.lazy_process and preprocessed:
            self.Warning.double_baseline()
        elif not self.lazy_process and not preprocessed:
            self.Warning.no_preprocessing()
        if self.titration is None:
            self.Warning.no_titration()
            return None
        try:
            spectra_wavelengths(self.data)
            names = spectrum_names(self.data)
        except InvalidWavelength:
            self.Error.invalid_wavelength()
            return None
        except SpectraError:
            self.Error.missing_spectrum_meta()
            return None

        rows = {name: index for index, name in enumerate(names)}
        buffer_row = rows.get(f"{BACKGROUND_PREFIX}buffer")
        sol_a_row = rows.get(f"{BACKGROUND_PREFIX}sol_A")
        sol_b_row = rows.get(f"{BACKGROUND_PREFIX}sol_B")
        for role, row in (("buffer", buffer_row), ("sol_A", sol_a_row)):
            if row is None:
                self.Error.missing_background(f"{BACKGROUND_PREFIX}{role}")
                return None
        if sol_b_row is None:
            self.Warning.no_sol_b()

        data_rows = [
            index
            for index, name in enumerate(names)
            if not name.startswith(BACKGROUND_PREFIX)
            and (split_series_name(name) or (None, None))[1] == DATA_STAGE
        ]
        if not data_rows:
            self.Error.no_data_spectra()
            return None

        try:
            ratio = self._titration_column(RATIO_COLUMN)
            if len(ratio) != len(data_rows):
                raise ValueError(
                    f"The Titration Table has {len(ratio)} rows but there "
                    f"are {len(data_rows)} '{DATA_STAGE}' spectra"
                )
            if DILUTION_COLUMN in self.titration.domain:
                dilution = self._titration_column(DILUTION_COLUMN)
            else:
                # Fixed-volume titrations do not dilute the cell.
                self.Warning.no_dilution_factor()
                dilution = np.ones(len(ratio))
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
        if self.lazy_process:
            X, _, (start, stop) = auto_baseline(X, self.lazy_window_percent)
            axis = spectra_wavelengths(self.data)[start:stop]
            region = (float(axis.min()), float(axis.max()))
            symbol = unit_symbol(table_unit(self.data, WAVELENGTH_UNIT_KEY)) or ""
            self.lazy_label.setText(
                f"Baseline taken from the flat region {region[0]:g}\u2013"
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

        samples = [split_series_name(names[row])[0] for row in data_rows]
        new_names, new_rows = [], []
        for stage in ("sol_A_buffer_subtracted", "sol_B_buffer_subtracted"):
            if stage in corrected:
                new_names.append(f"{BACKGROUND_PREFIX}{stage}")
                new_rows.append(corrected[stage][None, :])
        for stage in DATA_STAGES:
            new_names.extend(f"{sample} | {stage}" for sample in samples)
            new_rows.append(corrected[stage])

        # Keep the input rows so the original stages remain selectable.
        output = build_spectra_table(
            self.data,
            np.vstack((self.data.X, *new_rows)),
            [*names, *new_names],
            [*row_units, *[unit] * len(new_names)],
        )
        if region is not None:
            output.attributes["lazy_baseline_region"] = list(region)
        output.name = (
            f"{self.data.name} - corrected" if self.data.name else "Corrected CD data"
        )
        return output


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWCDDataCorrection).run()
