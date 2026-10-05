# CD Spectra Plot

The CD Spectra plot widget can be used to visualise sets of raw or processed spectra. 

It has been principally designed for the output of the [CD Data Loader](cd_data_loader.md) widget, or of spectral preprocessing applied to it.

Spectra are named `sample | stage`. The stages present in the table are read from the table and offered in the "Processing stage" selector, so subsets of the data can be readily selected and plotted. For example, selecting `raw_data` displays all of the titration data files, and `sol_A` the Solution A background spectrum.

Alternatively, groups of files can be interactively selected from the list provided in the "Spectra" window. 

The colour scale and line widths used in the plots can be interactively changed in the widget.

The widget does not generate any output, but the plot can be saved in the usual way for Orange.
