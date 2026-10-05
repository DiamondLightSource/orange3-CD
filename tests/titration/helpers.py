"""Shared builders for the widget tests."""

import numpy as np
from Orange.data import ContinuousVariable, Domain, StringVariable, Table

# Descending, as an instrument writes them.
WAVELENGTHS = [260.0, 250.0, 240.0, 230.0]


def write_csv(path, values, wavelengths=WAVELENGTHS):
    """Write a CD file in the layout read by ``file_parser``."""
    lines = ["Data:", "CircularDichroism,", "Wavelength,CD 1,CD 2", "nm,mdeg,mdeg"]
    lines += [f"{w},{v},{v}" for w, v in zip(wavelengths, values)]
    lines += ["", "HT,", "Wavelength,HT", "nm,V", "260,1", "250,1", "", ""]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def spectra_table(
    names,
    X,
    units=None,
    wavelengths=(230.0, 240.0, 250.0, 260.0),
    spectrum_unit="millidegree",
    wavelength_unit="nanometer",
):
    """A spectra-as-rows table; ``units`` adds a per-row ``Unit`` meta."""
    metas = [StringVariable("Spectrum")]
    columns = [np.asarray(names, dtype=object)]
    if units is not None:
        metas.append(StringVariable("Unit"))
        columns.append(np.asarray(units, dtype=object))
    table = Table.from_numpy(
        Domain([ContinuousVariable(str(w)) for w in wavelengths], metas=metas),
        np.asarray(X, dtype=float),
        metas=np.column_stack(columns),
    )
    if spectrum_unit:
        table.attributes["spectrum_unit"] = spectrum_unit
    if wavelength_unit:
        table.attributes["wavelength_unit"] = wavelength_unit
    return table


def titration_table(ratios, concentration=15.0, dilution=None):
    """Titration calculator style output with concentration unit on a column."""
    columns = [("normalised_molar_ratio", np.asarray(ratios, dtype=float))]
    if dilution is not None:
        columns.insert(0, ("dilution_factor", np.asarray(dilution, dtype=float)))
    columns.append(("working_concentration_a", np.full(len(ratios), concentration)))
    domain = Domain([ContinuousVariable(name) for name, _ in columns])
    domain["working_concentration_a"].attributes["unit"] = "micromolar"
    return Table.from_numpy(domain, np.column_stack([c for _, c in columns]))
