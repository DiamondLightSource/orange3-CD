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

|Variable|Description|
|-|-|
|`Ratio`| The input Molar Ratios in the Titration Calculator widget |
|`Predicted Volume`| The volume to add at this stage of the titration|
|`Volume Added This Step`| The volume to add at this stage of the titration|
|`Total Stock B Volume`| The total volume of Stock B that has been added to the cell at each titration point. |
|`Total Cell Volume`| The total volume of solution in the cell at each titration point. |
|`Dilution Factor`| The factor by which the initial cell volume has been diluted through the addition of the total volume of Stock B|
|`Normalised Molar Ratio`| The molar ratio of each titration point divided by the molar equivilent  of Stock B.|
|`Mode`| Will contain the same value on each row. Indicates whether the table has been set up for an experiment either in `Fixed` or `Increasing` mode. |
|`Volume Solution A`| Will contain the same value on each row. The initial volume of Solution A that is added to the cell at the start of the experiment. |
|`Volume Buffer`| Will contain the same value on each row. The initial volume of the buffer solution that is added to the cell. |
|`Max Volume Allowed`| Will contain the same value on each row. The maximum volume that can be added to the cell|
|`Max Volume Added`| Will contain the same value on each row. The total volume of Solution B that has been added by the end of the titration. |
|`Within Limit`| Should be `True` or `False`. Indicates whether the total volume of solution B added is less than the maximum volume of the cell. |

