#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Apply the titration reference corrections to loaded CD spectra."""

from __future__ import annotations

import numpy as np
from Orange.data import Table
from Orange.widgets import gui
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from . import Q_
from .utils import (
    BACKGROUND_PREFIX,
    InvalidWavelength,
    SpectraError,
    build_spectra_table,
    spectrum_names,
    spectrum_units,
    split_series_name,
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
    result.update({
        "buffer_subtraction": buffer_subtraction,
        "sol_A_subtraction": sol_a_subtraction,
        "subtract_frac_sol_B": subtract_frac_sol_b,
        "plus_sol_A": subtract_frac_sol_b + sol_a_bs,
    })
    return result


class OWCDDataCorrection(OWWidget):
    name = "CD Data Correction"
    description = (
        "Apply buffer, Solution A and Solution B corrections to preprocessed "
        "CD titration spectra."
    )
    icon = "icons/Titration.svg"
    priority = 25
    want_main_area = False
    resizing_enabled = False

    class Inputs:
        data = Input("CD Data", Table)
        titration = Input("Titration Table", Table)

    class Outputs:
        data = Output("Corrected CD Data", Table)

    class Error(OWWidget.Error):
        missing_spectrum_meta = Msg(
            "Input does not contain a 'Spectrum' string meta naming the rows."
        )
        invalid_wavelength = Msg(
            "Input attribute names are not all valid wavelengths."
        )
        missing_background = Msg("Background spectrum '{}' was not found.")
        no_data_spectra = Msg(f"No '{DATA_STAGE}' spectra were found.")
        invalid_titration = Msg("{}")
        incompatible_units = Msg("{}")

    class Warning(OWWidget.Warning):
        no_titration = Msg("Connect a Titration Table to correct the spectra.")
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
            "1. buffer subtraction x dilution factor\n"
            "2. subtract Solution A\n"
            "3. subtract Solution B x molar ratio\n"
            "4. add Solution A back\n"
            "The zero-level offset is not applied; use spectral preprocessing "
            "first.",
        )
        note.setWordWrap(True)

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
        if self.data is None:
            return None
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
            index for index, name in enumerate(names)
            if not name.startswith(BACKGROUND_PREFIX)
            and (split_series_name(name) or (None, None))[1] == DATA_STAGE
        ]
        if not data_rows:
            self.Error.no_data_spectra()
            return None

        try:
            dilution = self._titration_column(DILUTION_COLUMN)
            ratio = self._titration_column(RATIO_COLUMN)
            if len(dilution) != len(data_rows):
                raise ValueError(
                    f"The Titration Table has {len(dilution)} rows but there "
                    f"are {len(data_rows)} '{DATA_STAGE}' spectra"
                )
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

        def values(row: int) -> np.ndarray:
            return Q_(self.data.X[row], row_units[row]).to(unit).magnitude

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
        output.name = (
            f"{self.data.name} - corrected" if self.data.name else "Corrected CD data"
        )
        return output


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWCDDataCorrection).run()
