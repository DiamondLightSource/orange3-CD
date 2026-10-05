# CD Titration Calculator

Widget to generate titration input table.

On the boxes on the left of the widget in the `Titration inputs` section, add the following information about your experiment:
 - Stock solution concentrations and volumes.
 - Pipetting target range:
    * by setting minimum/maximum values, the calculator will select the best solution of stock B to add for your range.
 - Molar target ratios at which to titrate
    * this is a comma separated list of values, ie. `0.1, 0.2, 0.3...`
 - The titration mode.
    * This is either `Fixed` or `Increasing`, as selected from the drop down menu.


Once all this information is added, you can click the `Calculate` button on the bottom left of the widget, which will generate a table. 

The table is previewed on the right hand side of the widget. It may also be instructive to feed the table to the `Data Table` widget, where it can be previewed in full.

## Output table columns

Only quantities that differ between titration points are columns. Units are stored on the columns.

|Variable|Mode|Description|
|-|-|-|
|`ratio`|both| The input Molar Ratios in the Titration Calculator widget |
|`stock_b`|both| The number of the Stock B solution used at this point |
|`volume_added_this_step`|increasing| The volume of Stock B to add at this stage of the titration|
|`total_stock_b_volume`|increasing| The total volume of Stock B that has been added to the cell at each titration point. |
|`total_cell_volume`|increasing| The total volume of solution in the cell at each titration point. |
|`dilution_factor`|increasing| The factor by which the initial cell volume has been diluted through the addition of the total volume of Stock B|
|`volume_stock_b`|fixed| The volume of Stock B to add at this titration point|
|`baseline_volume`|fixed| The volume of buffer needed to bring the cell to its fixed volume|
|`concentration_b`|fixed| The concentration of Solution B in the cell at this point|
|`normalised_molar_ratio`|both| The molar ratio of each titration point divided by the molar equivalent of Stock B.|
|`working_concentration_a`|both| The working concentration of Solution A, read by the Delta Epsilon widget. |

Values that are the same for every titration point are stored as attributes of the table, not as columns:

|Attribute|Mode|Description|
|-|-|-|
|`mode`|both| Whether the table has been set up for an experiment in `fixed` or `increasing` mode. |
|`volume_solution_a`|both| The initial volume of Solution A that is added to the cell at the start of the experiment. |
|`volume_buffer`|increasing| The initial volume of the buffer solution that is added to the cell. |
|`max_volume_allowed`|increasing| The maximum total volume of Stock B that can be added to the cell (15% of the starting cell volume). |
|`max_volume_added`|increasing| The total volume of Stock B that has been added by the end of the titration. |
|`within_limit`|increasing| Whether the total volume of Stock B added is within the maximum allowed. A warning is shown in the widget when it is not. |
