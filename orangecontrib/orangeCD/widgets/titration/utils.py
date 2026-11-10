"""Helpers for the spectra-as-rows table layout shared by the CD widgets.

A spectra table has one row per spectrum and one continuous
attribute per wavelength, the attribute name being the wavelength itself.
Each row is labelled by string metas:

* ``Spectrum``: the name shown to the user, ``sample | stage``;
* ``Sample``: the file the spectrum came from, or the name of the reference
  (``buffer``, ``sol_A``, ``sol_B``);
* ``Stage``: how far the spectrum has been processed (``raw_data``,
  ``buffer_subtracted`` and so on);
* ``Role``: what the spectrum is, one of ``ROLES``.

Widgets select spectra by ``Sample``, ``Stage`` and ``Role`` and never take the
``Spectrum`` name apart. Tables with only a ``Spectrum`` meta, for example from
other add-ons, are still read by splitting the name.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from pint import Quantity

from ...units import Q_

SPECTRUM_META = "Spectrum"
UNIT_META = "Unit"
SAMPLE_META = "Sample"
STAGE_META = "Stage"
ROLE_META = "Role"
# Number of the titration point (1 for the first row of the Titration Table)
# that a data spectrum belongs to. NaN for spectra that are not titration points.
POINT_META = "Point"

# What a spectrum is. The three references are measured once; every other
# spectrum is a titration point.
ROLE_BUFFER = "buffer"
ROLE_SOL_A = "sol_A"
ROLE_SOL_B = "sol_B"
ROLE_TITRATION = "titration"
REFERENCE_ROLES = (ROLE_BUFFER, ROLE_SOL_A, ROLE_SOL_B)
ROLES = (*REFERENCE_ROLES, ROLE_TITRATION)

# How far a spectrum has been processed. Each stage is named for the step that
# produced it, and the same name is used for the references and the data.
STAGE_RAW = "raw_data"
STAGE_BUFFER_SUBTRACTED = "buffer_subtracted"
STAGE_SOL_A_SUBTRACTED = "sol_A_subtracted"
STAGE_SOL_B_SUBTRACTED = "sol_B_subtracted"
STAGE_SOL_A_ADDED = "sol_A_added"
DELTA_EPSILON_SUFFIX = "delta_epsilon"

# Table attribute keys. Units are stored as pint-parseable strings.
WAVELENGTH_UNIT_KEY = "wavelength_unit"
SPECTRUM_UNIT_KEY = "spectrum_unit"
CONCENTRATION_KEY = "solution_a_concentration"
PATHLENGTH_KEY = "pathlength"
SOLUTION_A_MW_KEY = "solution_a_molecular_weight"
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


def has_spectroscopy_preprocessing(table: Table) -> bool:
    """Whether Quasar (orangecontrib.spectroscopy) preprocessing was applied.

    Orange does not tell a widget which widget its input came from, but
    preprocessed attributes remember how they were computed. Quasar's
    preprocessors leave a ``compute_value`` from the
    ``orangecontrib.spectroscopy`` package on the attributes. A purely
    selecting step such as *Cut* leaves nothing behind, so it is not
    detected, and neither is a table that was saved and reloaded.
    """
    for variable in table.domain.attributes:
        compute = variable.compute_value
        for _ in range(50):  # guard against cyclic chains
            if compute is None:
                break
            if type(compute).__module__.startswith("orangecontrib.spectroscopy"):
                return True
            compute = getattr(getattr(compute, "variable", None), "compute_value", None)
    return False


def spectrum_names(table: Table) -> list[str]:
    """Return the name of every row, in row order."""
    for index, variable in enumerate(table.domain.metas):
        if variable.name == SPECTRUM_META and isinstance(variable, StringVariable):
            return [str(value) for value in table.metas[:, index]]
    raise MissingWavelength(SPECTRUM_META)


def wavelengths(table: Table) -> np.ndarray:
    """Return the wavelength axis, read from the attribute names."""
    spectrum_names(table)  # validates the layout
    try:
        return np.array([float(variable.name) for variable in table.domain.attributes])
    except ValueError as exc:
        raise InvalidWavelength(str(exc)) from exc


@dataclass(frozen=True)
class SpectrumLabel:
    """What one row of a spectra table is."""

    sample: str
    stage: str
    role: str

    @property
    def name(self) -> str:
        """The ``Spectrum`` name."""
        return f"{self.sample} | {self.stage}"


def _legacy_label(name: str) -> SpectrumLabel:
    """Read the role and stage from a ``Spectrum`` name alone.

    ``Background | sol_A`` is the Solution A reference, and
    ``Background | sol_A_buffer_subtracted`` the same after buffer subtraction.
    Anything else with a stage is a titration point.
    """
    parsed = split_series_name(name)
    if parsed is None:
        return SpectrumLabel(name, "", ROLE_TITRATION)
    sample, stage = parsed
    if sample == "Background":
        for role in REFERENCE_ROLES:
            if stage == role:
                return SpectrumLabel(role, STAGE_RAW, role)
            if stage.startswith(f"{role}_"):
                return SpectrumLabel(role, stage[len(role) + 1 :], role)
    return SpectrumLabel(sample, stage, ROLE_TITRATION)


def _string_meta(table: Table, name: str) -> list[str] | None:
    for index, variable in enumerate(table.domain.metas):
        if variable.name == name and isinstance(variable, StringVariable):
            return [str(value) for value in table.metas[:, index]]
    return None


def spectrum_labels(table: Table) -> list[SpectrumLabel]:
    """Return the sample, stage and role of every row, in row order."""
    names = spectrum_names(table)  # validates the layout
    samples = _string_meta(table, SAMPLE_META)
    stage_names = _string_meta(table, STAGE_META)
    if samples is None or stage_names is None:
        return [_legacy_label(name) for name in names]
    roles = _string_meta(table, ROLE_META) or [ROLE_TITRATION] * len(names)
    return [
        SpectrumLabel(sample, stage, role)
        for sample, stage, role in zip(samples, stage_names, roles)
    ]


def matching_spectra(
    table: Table, stage: str | None = None, role: str | None = None
) -> list[tuple[int, str]]:
    """Return ``(row, name)`` pairs for spectra of a stage and role.

    With ``stage=None`` every stage is returned, and with ``role=None`` every
    role.
    """
    names = spectrum_names(table)
    return [
        (row, names[row])
        for row, label in enumerate(spectrum_labels(table))
        if (stage is None or label.stage == stage)
        and (role is None or label.role == role)
    ]


def stages(table: Table, role: str | None = None) -> list[str]:
    """Return the stages present, in order of appearance."""
    found: list[str] = []
    for label in spectrum_labels(table):
        if label.stage and label.stage not in found and role in (None, label.role):
            found.append(label.stage)
    return found


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


def spectrum_points(table: Table) -> np.ndarray | None:
    """Return the titration point number of every row, or None if unlabelled.

    Rows that are not titration points have NaN. A table without the
    ``Point`` meta (for example one that was not made with a Titration Table
    connected to the CD Data Loader) returns None.
    """
    for index, variable in enumerate(table.domain.metas):
        if variable.name == POINT_META and isinstance(variable, ContinuousVariable):
            return np.asarray(table.metas[:, index], dtype=float)
    return None


def titration_rows(points: np.ndarray, count: int, what: str) -> np.ndarray:
    """Convert titration point numbers to zero-based Titration Table rows.

    ``points`` are the numbers of the spectra being paired, ``count`` the
    number of rows of the Titration Table. Raises ``ValueError`` if a spectrum
    has no point, a point does not exist, or two spectra share a point.
    """
    if np.isnan(points).any():
        raise ValueError(f"Some {what} spectra have no titration point")
    rows = points.astype(int) - 1
    if np.any(rows < 0) or np.any(rows >= count):
        raise ValueError(
            f"The {what} spectra use titration points that are not among "
            f"the {count} rows of the Titration Table"
        )
    if len(set(rows.tolist())) != len(rows):
        raise ValueError(f"Two {what} spectra have the same titration point")
    return rows


def build_spectra_table(
    template: Table,
    X: np.ndarray,
    names: list[SpectrumLabel] | list[str],
    units: list[str],
    attributes: dict | None = None,
    points: np.ndarray | None = None,
) -> Table:
    """Create a spectra table on the wavelength axis of ``template``.

    ``X`` holds one spectrum per row. ``names`` are either a ``SpectrumLabel``
    for each row, which give the ``Spectrum``, ``Sample``, ``Stage`` and
    ``Role`` metas, or plain names, which give only ``Spectrum``. ``units``
    gives the ``Unit`` meta. ``points``, if given, labels each row with its
    titration point (NaN for none) through the ``Point`` meta.
    """
    spectrum_variable = next(
        variable for variable in template.domain.metas if variable.name == SPECTRUM_META
    )
    labelled = bool(names) and isinstance(names[0], SpectrumLabel)
    meta_variables = [spectrum_variable]
    meta_columns = [
        np.asarray([n.name if labelled else n for n in names], dtype=object)
    ]
    if labelled:
        for meta, attribute in (
            (SAMPLE_META, "sample"),
            (STAGE_META, "stage"),
            (ROLE_META, "role"),
        ):
            meta_variables.append(StringVariable(meta))
            meta_columns.append(
                np.asarray([getattr(n, attribute) for n in names], dtype=object)
            )
    meta_variables.append(StringVariable(UNIT_META))
    meta_columns.append(np.asarray(units, dtype=object))
    if points is not None:
        meta_variables.append(ContinuousVariable(POINT_META))
        meta_columns.append(np.asarray(points, dtype=float).astype(object))
    domain = Domain(template.domain.attributes, metas=meta_variables)
    return Table.from_numpy(
        domain,
        np.asarray(X, dtype=float),
        metas=np.column_stack(meta_columns),
        attributes=dict(template.attributes if attributes is None else attributes),
    )
