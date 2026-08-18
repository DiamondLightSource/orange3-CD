# CD Titration Processing

Widget to process titration spectra.

This widget expects an input from the [CD Titration Calculator](titration_calculator.md) widget, so that the experimental details are used to process the spectra.


## Inputs
### Reference Files

The widget requires three reference spectra to process the data:
 - Solution A:
    * A spectrum of Solution A
 - Solution B:
    * A spectrum of Solution B
 - Buffer:
    * A spectrum of the buffer used in the experiment

The appropriate files can be selected by using the file browser from the `Browse...` button.

### Titration data files

The `Titration data files` section is used to load the data files from the experiment. The file browser is opened using the `Select files...` button.

``` 
NOTE

The loaded data files are sorted according to their name. 
```
 

### Background subtraction

A background subtraction region is defined by the text input boxes. The location of the region is indicated on the spectra on the right hand side of the widget by two vertical lines. The background region can also be set by moving the position of the lines with the mouse.


```
NOTE

The background subtraction is performed by taking a simple numeric average across the region set as the background for each spectra, and offsetting the entire spectrum by this value. This may change in future.
```


## Technical details

### Reference file background subtraction

The first stage of the processing is to process the reference files. The pipeline works as:

1) Subtract the buffer from both Solution A and Solution B.
    - These files are designated by `sol_A_buffer_subtracted` and `sol_B_buffer_subtracted` respectively
2) Make a zero offset to the buffer subtracted Solution A spectrum.
    - As described above, the background subtraction is performed by averaging over a certain region of the buffer subtracted spectrum, then subtracting this value from the entire spectrum.
    - This file is referred to as `sol_A_buffer_subtracted_zeroed`

### Data file processing

After processing the reference files, each data file is processed by the following pipeline:

1) Subtract the buffer from the spectrum
2) Scale the spectrum from 1) by the dilution factor of the titration point.
3) Subtract the `sol_A_buffer_subtracted` spectrum from above 
4) Take `sol_B_buffer_subtracted` and scale it by the normalised molar ratio of the titration point. Then subtract this spectrum from 3).
5) Make a zero offset to 4) 
6) Add `sol_A_buffer_subtracted_zeroed` to 5)


### Output

The processed output generates a data table that can best be viewed by the Data Table widget, or by the [CD Spectra Plot](./cd_spectra_plot.md) widget.

#### Data Table output
Viewing the output table in the Data Table widget, the columns are named according to both the input file name or designation, and the processing stage. In general:

|Column name| Description |
|-|-|
|Background files|
|`Background \| sol_A`| The raw spectrum of Solution A used as an input|
|`Background \| sol_B`| The raw spectrum of Solution B used as an input|
|`Background \| buffer`| The raw spectrum of the buffer used as an input|
|`Background \| sol_A_buffer_subtracted`| The buffer subtracted spectrum of solution A|
|`Background \| sol_B_buffer_subtracted`| The buffer subtracted spectrum of solution B|
|`Background \| sol_A_buffer_subtracted_zeroed`| The buffer subtracted spectrum of solution A, offset to zero|
|Data files. These are the same for all the files loaded as data, as indicated by the `*`. |
|`file * \| raw_data`| The raw spectrum of the data file used as an input|
|`file * \| buffer_subtraction`| The result of processing stage 2 above|
|`file * \| sol_A_subtraction`| The result of processing stage 3 above|
|`file * \| subtract_frac_sol_B`| The result of processing stage 4 above|
|`file * \| set_to_zero `| The result of processing stage 5 above|
|`file * \| plus_sol_A`| The result of processing stage 6 above|


#### CD Spectra Plot

The output can also be viewed in the [CD Spectra Plot](./cd_spectra_plot.md) widget.


