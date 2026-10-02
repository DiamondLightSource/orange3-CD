# CD Data Correction

Widget to apply the titration reference corrections to CD spectra. It is intended to be used **after** spectral preprocessing (for example baseline subtraction in Quasar's *Preprocess Spectra* widget):

`CD Data Loader` → `Preprocess Spectra` → `CD Data Correction` → [Delta Epsilon](delta_epsilon.md)

## Inputs

- **CD Data**: a spectra table from the [CD Data Loader](cd_data_loader.md), optionally preprocessed. It must contain the `Background | buffer` and `Background | sol_A` spectra, and the `raw_data` spectra of the titration points. `Background | sol_B` is optional; without it the Solution B subtraction is skipped and a warning is shown.
- **Titration Table**: from the [CD Titration Calculator](titration_calculator.md), providing the `dilution_factor` and `normalised_molar_ratio` of each titration point. It must have one row per `raw_data` spectrum.

All spectra are converted to the unit of the first data spectrum before the corrections, so spectra in different CD units (for example degree and millidegree) can be combined.

## Corrections

The reference spectra are first buffer subtracted:

- `Background | sol_A_buffer_subtracted` = Solution A − buffer
- `Background | sol_B_buffer_subtracted` = Solution B − buffer

Each data spectrum is then processed in turn:

1) `buffer_subtraction`: subtract the buffer, then scale by the dilution factor of the titration point.
2) `sol_A_subtraction`: subtract `sol_A_buffer_subtracted` from 1).
3) `subtract_frac_sol_B`: scale `sol_B_buffer_subtracted` by the normalised molar ratio of the titration point, and subtract it from 2).
4) `plus_sol_A`: add `sol_A_buffer_subtracted` back to 3).

> [!NOTE]
> The zero-level (baseline) offset that used to be part of the processing is not applied. Apply it with spectral preprocessing before this widget.

## Output

The output table contains all of the input rows, followed by the new rows named `Background | <stage>` and `<file> | <stage>`. Each row has a `Unit` meta. The stages can be selected in the [CD Spectra Plot](cd_spectra_plot.md) widget.

[Delta Epsilon](delta_epsilon.md) and Binding Data prefer the `plus_sol_A` and `sol_A_buffer_subtracted` series when they are present, and otherwise fall back to `raw_data` and `sol_A`.
