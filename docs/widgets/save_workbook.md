# Save workbook

This widget accepts an arbitrary number of input Data Tables as inputs. That is, any widget that generates tabulated data can be connected as an input.

The collection of inputs may then be saved as a single Excel workbook, with each input Data Table as a separate sheet. 

The save destination of the output data table is determined by the 'Save Workbook' button on the widget.

The names of the sheets can be interactively edited in the widget pane.


Column headers include the unit of the column where one is known, for example `Conc [B] (M)`. The columns of a spectra table are named by wavelength, so they carry the wavelength unit, for example `250.0 (nm)`.

If any of the input tables carry attributes, such as the Solution A concentration, the pathlength, the measurement wavelength or the fit settings, they are written to an extra `Metadata` sheet with the columns `Sheet`, `Attribute` and `Value`. The sheet is left out when there is nothing to record.

If there are no tables to save, or the workbook cannot be written (for example because the folder does not exist or the file is open in Excel), the widget shows an error instead of failing.
