# Secondary Structure

Estimates protein secondary structure fractions from delta epsilon CD spectra, using the CDPro methods (via the `cdpro` Python library).

## Input

A data table with a wavelength column (a continuous variable named `Wavelength`) and one or more delta epsilon columns. Delta epsilon columns are found by name (containing `delta_epsilon`) or by carrying the delta epsilon unit, so the output of the [Delta Epsilon](delta_epsilon.md) widget can be connected directly. Each delta epsilon column is fitted as a separate spectrum.

The data must have a value at every whole nanometre in its range. Rows with missing values are ignored for that spectrum, and non-integer wavelengths (e.g. 0.5 nm steps) are thinned to whole nanometres.

## Controls

- **Fitting method**: SELCON3, CDSSTR or CONTINLL.
- **Reference set**: the set of reference proteins, or Automatic, which picks one from the wavelength range.
- **CDSSTR random seed**: CDSSTR uses random subsets, so the seed is fixed to give repeatable results.
- **Wavelength limits**: optionally restrict the fit to a wavelength range.
- **Fit**: runs the fit in the background. Tick *Fit automatically* to refit on every change.

## Outputs

- **Secondary Structure**: one row per spectrum, with the structure fractions (names depend on the reference set), total helix and strand, RMSD and NRMSD of the fit. The Status and Message meta columns report spectra that could not be fitted.
- **Fitted Spectra**: the measured and calculated delta epsilon spectra, with a `Wavelength` meta column, suitable for plotting.
