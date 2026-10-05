# Secondary Structure

Estimates protein secondary structure fractions from delta epsilon CD spectra, using the CDPro methods (via the `cdpro` Python library).

## Input

A spectra table in the layout used by the other widgets: one row per spectrum, one attribute per wavelength, and a `Spectrum` string meta naming each row. Rows in delta epsilon units (from the `Unit` meta or `spectrum_unit` attribute) or with names ending `_delta_epsilon` are fitted, so the output of the [Delta Epsilon](delta_epsilon.md) widget can be connected directly. Wavelengths are converted to nanometres using the table's `wavelength_unit`.

The data must have a value at every whole nanometre in its range. Missing values are ignored for that spectrum, and non-integer wavelengths (e.g. 0.5 nm steps) are thinned to whole nanometres.

## Controls

- **Fitting method**: SELCON3, CDSSTR or CONTINLL.
- **Reference set**: the set of reference proteins, or Automatic, which picks one from the wavelength range.
- **CDSSTR random seed**: CDSSTR uses random subsets, so the seed is fixed to give repeatable results.
- **Wavelength limits**: optionally restrict the fit to a wavelength range.
- **Fit**: runs the fit in the background. Tick *Fit automatically* to refit on every change.

## Outputs

- **Secondary Structure**: one row per spectrum, with the structure fractions (names depend on the reference set), total helix and strand, RMSD and NRMSD of the fit. The Status and Message meta columns report spectra that could not be fitted.
- **Fitted Spectra**: a spectra table with `sample | measured` and `sample | calculated` rows in delta epsilon units, suitable for the CD Spectra Plot.
