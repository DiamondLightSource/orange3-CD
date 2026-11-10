# CD Data Loader

Widget to load circular dichroism spectra from CSV files into a table of spectra.

The widget applies **no corrections**. Its output is intended to feed the *Preprocess Spectra* widget from the Quasar spectroscopy add-on (`pip install orange3-cd[spectroscopy]`), where baseline subtraction and other preprocessing are done. The preprocessed table can then be passed to [Delta Epsilon](delta_epsilon.md) and the rest of the titration pipeline, or viewed with the [CD Spectra Plot](cd_spectra_plot.md) widget.

## Inputs

### Reference files (optional)

Spectra of Solution A, Solution B and the buffer can be loaded with the `Browse...` buttons. Any that are left empty are skipped. Solution A is required by [Delta Epsilon](delta_epsilon.md) and Binding Data downstream.

### Titration data files

The `Titration data files` section loads the data files from the experiment. The file browser is opened using the `Select files...` button.

> [!NOTE]
> The loaded data files are sorted according to their name. Every file in the list is always loaded. Selecting files in the list only picks which ones the `Move up` / `Move down` buttons move. The order of the files must match the order of the rows in the Titration Table.

### Titration Table (optional)

Connect the Titration Table from the [CD Titration Calculator](titration_calculator.md) to number the data files. The file at position *n* in the list is point *n* of the Titration Table. The list then reads, for example, `3. 0.3 equiv.csv  (ratio 0.3)`, so you can check at a glance that each file is paired with the right ratio, and reorder with `Move up` / `Move down` until it is.

Each data spectrum in the output then carries its point number in the `Point` meta. The later widgets use it to pick the right row of the Titration Table, instead of relying on the spectra still being in the same order as the table. Spectra that are not titration points (the references) have no point.

The number of files must equal the number of rows of the Titration Table. If it does not, a warning is shown and the files are not numbered. The later widgets then fall back to pairing by order.

### Units of the data files

The wavelength unit (nanometer, angstrom or micrometer) and CD signal unit (millidegree or degree) of the files are chosen here. They are not read from the files. The choice is stored on the output table, and every downstream widget uses it.

### Preview

The plot on the right shows the loaded data file spectra.

## Output

The output follows the layout used by the Quasar spectroscopy widgets: **each row is one spectrum and each column is one wavelength**. The column names are the wavelengths, in ascending order. The `Spectrum` meta column names each row:

|Spectrum| Description |
|-|-|
|`sol_A \| raw_data`| The spectrum of Solution A|
|`sol_B \| raw_data`| The spectrum of Solution B|
|`buffer \| raw_data`| The spectrum of the buffer|
|`file \| raw_data`| The spectrum of each titration data file, one row per file|

The name is for display only. What each row is, is recorded in separate meta columns, so a file name that happens to contain ` | ` is not misread:

|Meta|Description|
|-|-|
|`Sample`| The file name, or `sol_A`, `sol_B` or `buffer` for the references |
|`Stage`| How far the spectrum has been processed. The Loader writes `raw_data` |
|`Role`| What the spectrum is: `sol_A`, `sol_B`, `buffer` or `titration` |

If a Titration Table is connected there is also a `Point` meta, the number of the titration point of each data spectrum.

The units are stored on the table as the `wavelength_unit` and `spectrum_unit` attributes. If the files cover different wavelengths, only the wavelengths present in every file are kept and the widget shows a warning. The data is still sent. If the files have no wavelengths in common an error is shown instead.

> [!NOTE]
> Buffer subtraction, dilution-factor scaling and Solution B subtraction are not performed by this widget. Use the [CD Data Correction](cd_data_correction.md) widget after spectral preprocessing.
