# Glossary of terms, symbols and acronyms

This page explains the quantities, symbols and acronyms used by the widgets and in their documentation. The first section explains the physics, so that the names of the columns and spectra make sense. The sections after it are lookup tables.

## The physics

### Circular dichroism (CD)

Light can be circularly polarised, with its electric field rotating either clockwise or anticlockwise as it travels. Chiral molecules such as proteins absorb the two kinds of light by different amounts. Circular dichroism is that difference. A CD spectrum records it as a function of wavelength, and its shape reports on the structure of the molecule.

### Delta A (ΔA)

The difference in absorbance between left-handed and right-handed circularly polarised light:

ΔA = A<sub>L</sub> − A<sub>R</sub>

It has no unit. It is the quantity that the CD signal really measures, and it is the `Delta A` column of the [Binding Data](binding_data.md) table.

Like any absorbance, ΔA grows with the concentration of the sample and with the pathlength of the cell (see [Beer-Lambert law](#beer-lambert-law)).

### Millidegrees (mdeg) and 32980

CD instruments do not report ΔA. They report the *ellipticity* θ, an angle, usually in millidegrees (mdeg). The two are related by a fixed factor:

θ (mdeg) = 32980 × ΔA

so ΔA = θ / 32980. The number 32980 is 1000 × (ln 10 / 4) × (180 / π), where the 1000 is the conversion from degrees to millidegrees. It is the same for every sample and instrument. The widgets use it in [Delta Epsilon](delta_epsilon.md) and [Binding Data](binding_data.md) to turn a CD signal in millidegrees (or degrees) into ΔA.

### Beer-Lambert law

Absorbance is proportional to the concentration and the pathlength:

A = ε · c · l

where c is the concentration, l is the pathlength and ε is the molar absorption coefficient, a property of the molecule at that wavelength. The same law applies to the CD difference:

ΔA = Δε · c · l

### Delta epsilon (Δε)

The *molar circular dichroism coefficient*: how strongly one mole per litre of the molecule discriminates between left and right circularly polarised light, per centimetre of pathlength. Rearranging the law above:

Δε = ΔA / (c · l)

Its unit is L · mol<sup>-1</sup> · cm<sup>-1</sup>, which is the same as M<sup>-1</sup> cm<sup>-1</sup>.

The point of Δε is that it removes the dependence on concentration and pathlength. Two spectra recorded with different amounts of sample or in different cells cannot be compared as ΔA or millidegrees, but they can be compared as Δε.

### Mean-residue delta epsilon

For a protein, the secondary structure comes from the amide bonds of the chain, so what matters is the concentration of *residues* (amino acids), not of whole molecules. A protein of molecular weight MW<sub>A</sub> made of residues of average mass MRW contains MW<sub>A</sub> / MRW residues, so the residue concentration is c × MW<sub>A</sub> / MRW. The mean-residue delta epsilon is Δε calculated on that basis:

Δε<sub>MRW</sub> = ΔA × (MRW / MW<sub>A</sub>) / (c · l)

This is the quantity the [Delta Epsilon](delta_epsilon.md) widget produces and the one that secondary structure estimation uses, because it does not depend on the size of the protein. The typical value of the mean residue molecular weight (MRW) is about 113 g/mol, which is the widget's default.

The `Delta Epsilon` column of the [Binding Data](binding_data.md) table is different. It is the ordinary, per-molecule Δε = ΔA / (c · l), so the two are related by the factor MRW / MW<sub>A</sub>.

## The titration experiment

In a CD titration a solution of one molecule (the host) is held at a fixed concentration, and a second molecule (the ligand) is added to it in steps. The CD spectrum is recorded after each addition, and the change in the spectrum shows how much of the host has bound the ligand.

| Term | Meaning |
|-|-|
| **Solution A** | The molecule at a fixed concentration, normally the protein. Called *Host A* in CD Apps. |
| **Solution B** | The molecule that is added in steps, normally the ligand. Called *Ligand B* in CD Apps. |
| **Buffer** | The solution without A or B, measured as the blank. Called the *baseline* in CD Apps. It is subtracted from every spectrum. |
| **Stock solution** | A concentrated solution that is pipetted from. *Stock A* is concentrated Solution A, and *Stock B 1, 2, 3...* are Solution B at different concentrations, so that a convenient volume can be pipetted at each step. |
| **Working concentration of A** | The concentration of Solution A in the cell, as opposed to in its stock. |
| **Starting cell volume** | The volume in the cell before any Solution B is added. |
| **Titration point** | One step of the titration, which has one spectrum. |
| **Molar ratio** | The concentration of B divided by the concentration of A in the cell at a titration point. It is the `ratio` column of the Titration Table and the `Titration point` column of the Binding Data table. |
| **Stock B molar equivalent** | The molar ratio at which the reference spectrum of Solution B was recorded. |
| **Normalised molar ratio** | The molar ratio divided by the Stock B molar equivalent. It is the fraction of the Solution B reference spectrum that is present at a titration point, and is the `normalised_molar_ratio` column. |
| **Pathlength** | The distance the light travels through the sample, which is the width of the cell, in centimetres. |
| **Fixed volume** | Titration mode in which the cell always holds the same total volume. Buffer is removed to make room for each addition, so the concentration of A does not change. |
| **Increasing volume** | Titration mode in which each addition of B adds to the volume in the cell. A is progressively diluted, which has to be corrected for. |
| **Dilution factor** | In increasing mode, the total volume in the cell divided by the starting volume. It is 1 at the start, and slightly above 1 after each addition. Fixed volume titrations have no dilution. |
| **15% limit** | In increasing mode the total volume of B added should not exceed 15% of the starting cell volume, or the cell is diluted too much. The Titration Calculator warns if it is exceeded. Use a more concentrated Stock B or fewer or smaller steps. |
| **Equiv** | Short for *equivalents*. A file called `1.6 equiv.csv` is the titration point at a molar ratio of 1.6. |

## Spectra and their processing stages

The tables passed between the widgets have one row per spectrum. Each row has a `Sample` (the file, or the reference it is), a `Stage` (how far it has been processed) and a `Role` (what it is), and its display name is `sample | stage`. The same stage name means the same step whether it applies to a reference or to a titration spectrum.

**Roles**

| Role | Meaning |
|-|-|
| `buffer` | The buffer reference. |
| `sol_A` | The Solution A reference. |
| `sol_B` | The Solution B reference. |
| `titration` | A titration point. |

**Stages**

| Stage | Meaning |
|-|-|
| `raw_data` | The measured spectrum, averaged over its repeats. Every reference and every titration file starts here. |
| `buffer_subtracted` | The spectrum minus the buffer. For a titration spectrum it is also multiplied by the dilution factor. For a reference (`sol_A` or `sol_B`) it is simply the reference minus the buffer. |
| `sol_A_subtracted` | The `buffer_subtracted` titration spectrum minus the buffer-subtracted Solution A. This leaves only what the addition of B contributes. |
| `sol_B_subtracted` | `sol_A_subtracted` minus the buffer-subtracted Solution B scaled by the normalised molar ratio. This removes the part of the signal that is just the free ligand. What remains is the change caused by binding. |
| `sol_A_added` | `sol_B_subtracted` with the buffer-subtracted Solution A added back, so that the spectrum is that of A, changed by binding. |
| `<stage>_delta_epsilon` | The spectrum of that stage converted to mean-residue Δε by the Delta Epsilon widget. |

Other terms:

| Term | Meaning |
|-|-|
| **Average** | The mean of the repeated scans in a file, which the [CD Data Loader](cd_data_loader.md) reads as the spectrum. |
| **Zero-level offset** | A constant that is subtracted so that the spectrum is zero where the sample does not absorb, normally at long wavelengths. |
| **Lazy process** | The option in [CD Data Correction](cd_data_correction.md) that finds the flattest stretch of each spectrum and subtracts its mean as a zero-level offset. |
| **`Spectrum` meta** | The string column that names each row of a spectra table, `sample | stage`. It is for display: the widgets select spectra from the `Sample`, `Stage` and `Role` metas. Tables from other sources that have only a `Spectrum` meta are read by splitting the name. |
| **`Sample`, `Stage` and `Role` metas** | The string columns that say which file or reference a row is, how far it has been processed, and what it is. |
| **`Point` meta** | The number of the titration point (the row of the Titration Table) that a data spectrum belongs to. It is added by the CD Data Loader when the Titration Table is connected, and lets later widgets pair spectra with table rows by number rather than by order. |
| **`Unit` meta** | The string column that gives the unit of each row, for tables that mix millidegrees and Δε. |
| **Table attributes** | Values stored on a whole table rather than a column, such as `wavelength_unit`, `spectrum_unit`, `solution_a_concentration` and `pathlength`. They are how information travels from one widget to the next. |

## Columns of the Binding Data table

| Column | Meaning |
|-|-|
| `Titration point` | The molar ratio of B to A (0 for Solution A alone). |
| `CD` | The signal of the selected series at the chosen wavelength. |
| `Change in CD` | `CD` minus the CD of Solution A. By default the absolute value is taken. |
| `Delta A` | `Change in CD` expressed as ΔA, that is, divided by 32980 if it is in millidegrees. |
| `Delta Epsilon` | `Delta A / (concentration × pathlength)`. |
| `Binding Stoichiometry` | `ratio / (ratio + 1)`, the mole fraction of B in the mixture. Plotting the signal against it is a way to read off the binding stoichiometry, where a maximum at 0.5 suggests one A binds one B. |
| `Conc [B]` | The concentration of B in the cell, which is the molar ratio multiplied by the concentration of Solution A, in molar. This is what Origin plots. |
| `Hill Conc [B]` | `Conc [B]` in micromolar. |
| `Hill Delta A` | `Delta A` divided by the pathlength, which is ΔA per 1 cm. |

## Fitting

[Binding Plot](binding_plot.md) fits curves to the binding data. The curves describe how a signal y rises with x as more ligand is added and the binding sites fill.

| Symbol | Meaning |
|-|-|
| **Hill equation** | A curve that rises from zero to a maximum, `y = v_max · xⁿ / (K_half ⁿ + xⁿ)`. |
| `v_max` | The maximum of the curve, which is the signal when binding is complete (saturation). |
| `K_half` (`half_saturation`) | The value of x at which the signal is half of its maximum. It has the units of x, and a smaller value means tighter binding. |
| `n` (`hill_coefficient`) | The steepness of the curve. A value of 1 is simple one-site binding, above 1 suggests cooperative binding. |
| `bottom`, `top` | In the four-parameter Hill equation, the signal before any binding and at saturation. It replaces `v_max` and allows a starting offset. |
| **BiHill equation** | A curve that rises and then falls. |
| `p_m` | The amplitude of the BiHill curve. It is the height of the peak when the rise and the fall are well separated. |
| `k_a`, `h_a` | The position (half-way up) and steepness of the rising part of the BiHill curve. |
| `k_i`, `h_i` | The position and steepness of the falling part of the BiHill curve. |
| `r_sq` (R²) | The coefficient of determination. It is 1 for a perfect fit, and falls towards 0 as the fit gets worse. |
| `red_chi_sq` (reduced χ²) | The sum of squared residuals per degree of freedom. It is a measure of the scatter of the points around the curve, and is only comparable between fits of the same data. |
| **Standard error** | The uncertainty of a fitted parameter estimated by lmfit. |

## Units

| Symbol | Meaning |
|-|-|
| nm, Å, µm | Nanometre, ångström (0.1 nm) and micrometre, the wavelength units the Loader accepts. |
| mdeg, deg | Millidegree and degree, the units of the CD signal. |
| M, mM, µM, nM | Molar (mol/L), millimolar, micromolar and nanomolar. |
| L, mL, µL, nL | Litre, millilitre, microlitre and nanolitre. |
| g/mol | Grams per mole, the unit of molecular weight. |
| M⁻¹ cm⁻¹ | The unit of Δε. It is the same as L · mol⁻¹ · cm⁻¹. |

Units are handled by the `pint` library and are stored with the data, so every widget converts between them when needed.

## Software and file names

| Name | Meaning |
|-|-|
| **CD Apps** | The original Windows program from which these widgets are ported. |
| **Chirascan** | The CD spectrometer whose `.csv` files the [CD Data Loader](cd_data_loader.md) reads. |
| **Quasar** | The spectroscopy add-on for Orange (`orange-spectroscopy`). Its *Preprocess Spectra* widget can be used between the Loader and CD Data Correction, for example to subtract a baseline. |
| **Origin** | The data-analysis program that CD Apps exported its binding data to, which is why the Binding Data table has the columns it does. |
| **lmfit** | The Python library used to fit the Hill and BiHill equations. |
| **pint** | The Python library used for units. |
| **SSE** | Secondary structure estimation: working out how much of a protein is α-helix, β-sheet and other structure from its CD spectrum. See [Secondary Structure](secondary_structure.md). |
| **MRW** | Mean residue weight (molecular weight), the average mass of one amino acid residue, about 113 g/mol. |
| **MW<sub>A</sub>** | The molecular weight of Solution A. |
