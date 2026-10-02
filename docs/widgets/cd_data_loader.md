# CD Data Loader

Widget to load circular dichroism spectra from CSV files into a table of spectra.

The widget applies **no corrections**. Its output is intended to feed the *Preprocess Spectra* widget from the Quasar spectroscopy add-on (`pip install orange3-cd[spectroscopy]`), where baseline subtraction and other preprocessing are done. The preprocessed table can then be passed to [Delta Epsilon](delta_epsilon.md) and the rest of the titration pipeline, or viewed with the [CD Spectra Plot](cd_spectra_plot.md) widget.

## Inputs

### Reference files (optional)

Spectra of Solution A, Solution B and the buffer can be loaded with the `Browse...` buttons. Any that are left empty are skipped. Solution A is required by [Delta Epsilon](delta_epsilon.md) and Binding Data downstream.

### Titration data files

The `Titration data files` section loads the data files from the experiment. The file browser is opened using the `Select files...` button.

> [!NOTE]
> The loaded data files are sorted according to their name.

### Units of the data files

The wavelength unit (nanometer, angstrom or micrometer) and CD signal unit (millidegree or degree) of the files are chosen here. They are not read from the files. The choice is stored on the output table, and every downstream widget uses it.

### Preview

The plot on the right shows the loaded data file spectra.

## Output

The output follows the layout used by the Quasar spectroscopy widgets: **each row is one spectrum and each column is one wavelength**. The column names are the wavelengths, in ascending order. The `Spectrum` meta column names each row:

|Spectrum| Description |
|-|-|
|`Background \| sol_A`| The spectrum of Solution A|
|`Background \| sol_B`| The spectrum of Solution B|
|`Background \| buffer`| The spectrum of the buffer|
|`file \| raw_data`| The spectrum of each titration data file, one row per file|

The units are stored on the table as the `wavelength_unit` and `spectrum_unit` attributes. All files must share the same wavelength axis.

> [!NOTE]
> Buffer subtraction, dilution-factor scaling and Solution B subtraction are not performed by this widget. Use the [CD Data Correction](cd_data_correction.md) widget after spectral preprocessing.
