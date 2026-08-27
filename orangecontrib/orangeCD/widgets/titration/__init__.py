NAME = "CD Titration"

ID = "orangecontrib.orangeCD.widgets.titration"

DESCRIPTION = "Tools for CD titration data analysis."

ICON = "icons/Titration.svg"

BACKGROUND = "#C7D9F1"

PRIORITY = 100

# Use pint to handle units across the pipeline
from pint import UnitRegistry, set_application_registry
ureg = UnitRegistry()
Q_ = ureg.Quantity
set_application_registry(ureg)
