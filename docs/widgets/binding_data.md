# Binding Data

Extracts the CD signal of the titration at one wavelength and builds the binding table used for plotting and fitting (including the columns used for Origin / CD Apps).

## Inputs

- **CD and Delta Epsilon Spectra**: the output of the [Delta Epsilon](delta_epsilon.md) widget. The Solution A concentration and pathlength are read from the attributes it stores on the table.
- **Titration Table**: from the [CD Titration Calculator](titration_calculator.md). It must have one row per titration spectrum and a `ratio` (or `normalised_molar_ratio` / `molar_ratio`) column.

## Controls

- **Measurement wavelength**: the wavelength to extract. It can also be set by dragging the red line on the plot, and snaps to the nearest wavelength in the data.
- **Titration data series** and **Solution A series**: as in [Delta Epsilon](delta_epsilon.md). The series must be in CD units (for example millidegree), not delta epsilon.
- **Use absolute change from Solution A**: take the absolute value of the difference from Solution A. If unticked the signed difference is used.

The plot shows the Solution A and titration spectra and the selected wavelength.

## Output

**Binding Data**: one row for Solution A (titration point 0) followed by one row per titration point. The `Sample` meta gives the spectrum name.

|Column|Description|
|-|-|
|`Titration point`| The molar ratio (0 for Solution A) |
|`CD`| The CD signal at the selected wavelength |
|`Change in CD`| The difference from Solution A |
|`Delta A`| The change in CD as an absorbance difference |
|`Delta Epsilon`| `Delta A / (concentration × pathlength)` |
|`Binding Stoichiometry`| `titration point / (titration point + 1)` |
|`Conc [B]`| The titration point multiplied by the molar concentration of Solution A |

The table attributes record the measurement wavelength, the series used, the concentration and the pathlength. The output can be connected to [Binding Plot](binding_plot.md).
