"""Helpers for the spectra-as-rows table layout shared by the CD widgets.

A processed spectra table has one row per spectrum and one continuous
attribute per wavelength, the attribute name being the wavelength itself.
Each row is labelled by a string meta called ``Wavelength`` holding names of
the form ``sample | processing_stage``.
"""

from __future__ import annotations

import numpy as np
from Orange.data import StringVariable, Table

WAVELENGTH_META = "Wavelength"
BACKGROUND_PREFIX = "Background | "


class SpectraError(ValueError):
    """The table does not follow the expected spectra layout."""


class MissingWavelength(SpectraError):
    """The spectrum-name meta is absent or is not a string variable."""


class InvalidWavelength(SpectraError):
    """An attribute name could not be read as a wavelength."""


def split_series_name(name: str) -> tuple[str, str] | None:
    """Split a spectrum name of the form ``sample | processing_stage``."""
    if " | " not in name:
        return None
    sample, stage = name.rsplit(" | ", maxsplit=1)
    sample, stage = sample.strip(), stage.strip()
    return (sample, stage) if sample and stage else None


def spectrum_names(table: Table) -> list[str]:
    """Return the name of every row, in row order."""
    for index, variable in enumerate(table.domain.metas):
        if variable.name == WAVELENGTH_META and isinstance(
            variable, StringVariable
        ):
            return [str(value) for value in table.metas[:, index]]
    raise MissingWavelength(WAVELENGTH_META)


def wavelengths(table: Table) -> np.ndarray:
    """Return the wavelength axis, read from the attribute names."""
    spectrum_names(table)  # validates the layout
    try:
        return np.array(
            [float(variable.name) for variable in table.domain.attributes]
        )
    except ValueError as exc:
        raise InvalidWavelength(str(exc)) from exc


def matching_spectra(
    table: Table, stage: str | None = None
) -> list[tuple[int, str]]:
    """Return ``(row, name)`` pairs for spectra of a processing stage.

    With ``stage=None`` every spectrum is returned.
    """
    pairs = list(enumerate(spectrum_names(table)))
    if stage is None:
        return pairs
    suffix = f" | {stage}"
    return [(row, name) for row, name in pairs if name.endswith(suffix)]


def stages(table: Table) -> list[str]:
    """Return the processing stages present, in order of appearance."""
    found: list[str] = []
    for name in spectrum_names(table):
        parsed = split_series_name(name)
        if parsed and parsed[1] not in found:
            found.append(parsed[1])
    return found


def reference_spectrum(
    candidates: list[tuple[int, str]],
) -> tuple[int, str] | None:
    """Pick the Solution A reference: a Background spectrum, else the first."""
    for candidate in candidates:
        if candidate[1].startswith(BACKGROUND_PREFIX):
            return candidate
    return candidates[0] if candidates else None
