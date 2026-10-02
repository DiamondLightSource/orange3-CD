import os

# Widgets need a Qt platform; run headless unless the caller chose one.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest


@pytest.fixture(autouse=True)
def _expose_tmp_path(request, tmp_path):
    """Give unittest-style (WidgetTest) classes ``self.tmp_path``."""
    if request.instance is not None:
        request.instance.tmp_path = tmp_path
