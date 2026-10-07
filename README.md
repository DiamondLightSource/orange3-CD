# Orange3-CD

Orange3 add on for Circular Dichroism analysis. Currently under development.

The widgets are documented in [docs/widgets](docs/widgets). They are grouped in two categories: *CD Titration* (data loading, correction, titration calculation, binding data and fitting) and *CD Core* ([Secondary Structure](docs/widgets/secondary_structure.md) analysis).


### Current TODO

 - Finish documentation (Save Workbook and example pipelines)
 - Fitting:
   * Add additional fitting functional forms to the `Binding Plot` widget
   * Port the fitting to lmfit for nicer reporting

### Longer term TODO

 - Adopt the [DLS Python Copier Template](https://github.com/DiamondLightSource/python-copier-template/) for this project for full CD/CI
    * Because we have to still follow the guidelines for Orange3 packages, this must be done carefully
 - Convert the documentation to a readthedocs site
 - Draft example pipelines for different experimental set ups
 - Add more pipeline widgets

