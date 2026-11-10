# CD Data Correction

Widget to apply the titration reference corrections to CD spectra. It is intended to be used **after** spectral preprocessing (for example baseline subtraction in Quasar's *Preprocess Spectra* widget):

`CD Data Loader` → `Preprocess Spectra` → `CD Data Correction` → [Delta Epsilon](delta_epsilon.md)

## Inputs

- **CD Data**: a spectra table from the [CD Data Loader](cd_data_loader.md), optionally preprocessed. It must contain the raw `buffer` and `sol_A` reference spectra (the rows whose `Role` is `buffer` or `sol_A` and whose `Stage` is `raw_data`), and the `raw_data` spectra of the titration points (`Role` `titration`). The `sol_B` reference is optional; without it the Solution B subtraction is skipped and a warning is shown.
- **Titration Table**: from the [CD Titration Calculator](titration_calculator.md), providing the `dilution_factor` and `normalised_molar_ratio` of each titration point. If the spectra carry a `Point` meta (the Titration Table was connected to the [CD Data Loader](cd_data_loader.md)), each `raw_data` spectrum is paired with the table row of its point, whatever the order of the spectra, and a subset of the points can be used. Otherwise the table must have one row per `raw_data` spectrum, and they are paired in order. Fixed-volume tables have no `dilution_factor` column; the widget then assumes a factor of 1 for every point, shows a warning, and still sends the output.

All spectra are converted to the unit of the first data spectrum before the corrections, so spectra in different CD units (for example degree and millidegree) can be combined.

## Corrections

The reference spectra are first buffer subtracted:

- `sol_A | buffer_subtracted` = Solution A − buffer
- `sol_B | buffer_subtracted` = Solution B − buffer

Each data spectrum is then processed in turn:

1) `buffer_subtracted`: subtract the buffer, then scale by the dilution factor of the titration point.
2) `sol_A_subtracted`: subtract the buffer-subtracted Solution A from 1).
3) `sol_B_subtracted`: scale the buffer-subtracted Solution B by the normalised molar ratio of the titration point, and subtract it from 2).
4) `sol_A_added`: add the buffer-subtracted Solution A back to 3).

> [!NOTE]
> The zero-level (baseline) offset that used to be part of the processing is not applied. Apply it with spectral preprocessing before this widget.

## Output

The output table contains all of the input rows, followed by the new rows named `<sample> | <stage>`, with the `Sample`, `Stage` and `Role` metas set. The references are `sol_A` and `sol_B` with the stage `buffer_subtracted`. Each row has a `Unit` meta. The stages can be selected in the [CD Spectra Plot](cd_spectra_plot.md) widget.

[Delta Epsilon](delta_epsilon.md) and Binding Data prefer the `sol_A_added` titration series and the `buffer_subtracted` Solution A series when they are present, and otherwise fall back to `raw_data` for both.
