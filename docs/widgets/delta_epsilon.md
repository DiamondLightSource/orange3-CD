# Delta Epsilon

Converts CD spectra to mean-residue delta epsilon. The selected CD spectra are kept in the output, and a matching delta epsilon spectrum is added for each one.

`CD Data Correction` → `Delta Epsilon` → [Secondary Structure](secondary_structure.md) / Binding Data

## Inputs

- **CD Data**: a spectra table (one row per spectrum, one attribute per wavelength, a `Spectrum` meta naming each row as `sample | series`), normally the output of the [CD Data Correction](cd_data_correction.md) widget.
- **Titration Table**: from the [CD Titration Calculator](titration_calculator.md). Its `working_concentration_a` column provides the Solution A concentration. Without it a warning is shown and nothing is converted.

## Controls

- **Titration data series**: the series to convert, for example `plus_sol_A` or `raw_data`. Defaults to `plus_sol_A` when present, otherwise `raw_data`.
- **Solution A series**: the Solution A reference series to convert alongside it. Defaults to `sol_A_buffer_subtracted` when present, otherwise `sol_A`.
- **Pathlength** (cm), **Mean residue molecular weight** and **Solution A molecular weight** (g/mol): the conversion parameters. All must be greater than zero.
- **Apply**: sends the output. Tick *Apply automatically* to update on every change.

The conversion is

delta epsilon = (CD / 32980 mdeg) × mean residue molecular weight / (concentration × pathlength × Solution A molecular weight)

Input spectra in other CD units (for example degree) are converted to millidegree first.

## Output

**CD and Delta Epsilon Spectra**: the selected CD spectra, followed by rows named `<sample> | <series>_delta_epsilon` with their `Unit` meta set to delta epsilon. The concentration and pathlength used are stored on the table, so that [Binding Data](binding_data.md) can use them. The output can be connected to [Secondary Structure](secondary_structure.md), [Binding Data](binding_data.md) or the [CD Spectra Plot](cd_spectra_plot.md).
