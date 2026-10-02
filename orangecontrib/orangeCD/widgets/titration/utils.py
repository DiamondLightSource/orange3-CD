"""Helpers for the spectra-as-rows table layout shared by the CD widgets.

A processed spectra table has one row per spectrum and one continuous
attribute per wavelength, the attribute name being the wavelength itself.
Each row is labelled by a string meta called ``Spectrum`` holding names of
the form ``sample | processing_stage``.
"""

from __future__ import annotations

import numpy as np
from Orange.data import StringVariable, Table
from pint import Quantity

from . import Q_

SPECTRUM_META = "Spectrum"
UNIT_META = "Unit"
BACKGROUND_PREFIX = "Background | "

# Table attribute keys. Units are stored as pint-parseable strings.
WAVELENGTH_UNIT_KEY = "wavelength_unit"
SPECTRUM_UNIT_KEY = "spectrum_unit"
CONCENTRATION_KEY = "solution_a_concentration"
PATHLENGTH_KEY = "pathlength"
MEASUREMENT_WAVELENGTH_KEY = "measurement_wavelength"

# CD signal that corresponds to one unit of absorbance difference.
MDEG_PER_DELTA_A = Q_(32980.0, "millidegree")
DELTA_EPSILON_UNIT = "liter / mole / centimeter"


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
        if variable.name == SPECTRUM_META and isinstance(
            variable, StringVariable
        ):
            return [str(value) for value in table.metas[:, index]]
    raise MissingWavelength(SPECTRUM_META)


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


def unit_string(unit: str) -> str:
    """Normalise a unit name through pint, e.g. ``"nanometer"``."""
    return str(Q_(1, unit).units)


def unit_symbol(unit: str | None) -> str | None:
    """Short display symbol for a unit, e.g. ``"nm"``; None if no unit."""
    if not unit:
        return None
    return f"{Q_(1, unit).units:~}"


def table_unit(table: Table, key: str) -> str | None:
    """Return the unit string stored in ``table.attributes[key]``."""
    value = table.attributes.get(key)
    return str(value) if value else None


def quantity_string(quantity: Quantity) -> str:
    """Serialise a quantity for storage in ``table.attributes``."""
    return f"{float(quantity.magnitude)!r} {quantity.units}"


def table_quantity(table: Table, key: str) -> Quantity:
    """Read a quantity stored with :func:`quantity_string`.

    Raises ``KeyError`` if absent and ``ValueError`` if it cannot be parsed
    or is not a finite positive number.
    """
    text = table.attributes[key]
    try:
        quantity = Q_(text)
        magnitude = float(quantity.magnitude)
    except Exception as exc:  # pint raises several error types
        raise ValueError(f"{key}: cannot parse {text!r}") from exc
    if not np.isfinite(magnitude) or magnitude <= 0:
        raise ValueError(f"{key}: must be greater than zero")
    return quantity


def spectrum_units(table: Table) -> list[str | None]:
    """Return the unit of every spectrum row.

    A ``Unit`` string meta takes priority (needed when one table mixes
    units, e.g. mdeg and delta epsilon); otherwise the table-wide
    ``spectrum_unit`` attribute applies to every row.
    """
    for index, variable in enumerate(table.domain.metas):
        if variable.name == UNIT_META and isinstance(variable, StringVariable):
            return [str(value) or None for value in table.metas[:, index]]
    return [table_unit(table, SPECTRUM_UNIT_KEY)] * len(table)


def shared_unit(units) -> str | None:
    """Return the single unit common to all of ``units``, else None."""
    unique = set(units)
    return unique.pop() if len(unique) == 1 else None
