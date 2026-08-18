# CD Spectra Plot

The CD Spectra plot widget can be used to visualise sets of raw or processed spectra. 

It has been principally designed for the output of the [CD Titration Processing](titration_processing.md) widget. 

TODO: in the widget, the PROCESSING_STAGES variable is hard coded. This should really look at the table and determine what stages of processing are available first. Also needs to sensibly include the background spectrum for these data sets.

From the processing stage of the data as described in the documentation for the [CD Titration Processing](./titration_processing.md) widget, subsets of the data can be readily selected and plotted. 

For example, selecting the `plus_sol_A` processing stage will display the final processed spectra of all of the data.

Alternatively, groups of files can be interactively selected from the list provided in the "Spectra" window. 

The colour scale and line widths used in the plots can be interactively changed in the widget.

The widget does not generate any output, but the plot can be saved in the usual way for Orange.
