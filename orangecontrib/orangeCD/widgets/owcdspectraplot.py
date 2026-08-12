#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Orange widget for plotting processed circular-dichroism spectra."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from AnyQt.QtCore import Qt
from AnyQt.QtGui import QColor
from AnyQt.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QGroupBox,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QVBoxLayout,
)
from Orange.data import ContinuousVariable, Table
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, OWWidget

PROCESSING_STAGES = (
    "All spectra",
    "raw_data",
    "buffer_subtraction",
    "sol_A_subtraction",
    "subtract_frac_sol_B",
    "plus_sol_A",
)

COLOUR_SCALES = {
    "Viridis": ("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"),
    "Plasma": ("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921"),
    "Inferno": ("#000004", "#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4"),
    "Cividis": ("#00224e", "#31446b", "#666970", "#a38f63", "#e6c75a", "#fee838"),
    "Blue to red": ("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b"),
    "Greyscale": ("#111111", "#555555", "#999999", "#dddddd"),
    "Classic": ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"),
}


def colours_from_scale(scale_name: str, count: int) -> list[QColor]:
    """Sample *count* line colours across the selected colour scale."""
    if count <= 0:
        return []
    anchors = [QColor(value) for value in COLOUR_SCALES[scale_name]]
    if count == 1:
        return [anchors[len(anchors) // 2]]

    colours: list[QColor] = []
    for position in np.linspace(0.0, len(anchors) - 1, count):
        lower = int(np.floor(position))
        upper = min(lower + 1, len(anchors) - 1)
        fraction = position - lower
        first, second = anchors[lower], anchors[upper]
        colours.append(
            QColor.fromRgbF(
                first.redF() + fraction * (second.redF() - first.redF()),
                first.greenF() + fraction * (second.greenF() - first.greenF()),
                first.blueF() + fraction * (second.blueF() - first.blueF()),
            )
        )
    return colours


class OWCDSpectraPlot(OWWidget):
    name = "CD Spectra Plot"
    description = "Plot spectra produced by the CD Titration Processing widget."
    icon = "icons/Titration.svg"
    priority = 30
    want_main_area = True
    resizing_enabled = True

    selected_stage = Setting("plus_sol_A")
    colour_scale = Setting("Viridis")
    show_legend = Setting(True)
    line_width = Setting(1.5)
    reverse_wavelength_axis = Setting(False)

    class Inputs:
        data = Input("Processed CD Data", Table)

    class Error(OWWidget.Error):
        missing_wavelength = Msg(
            "The input table does not contain a continuous Wavelength meta attribute."
        )
        no_numeric_spectra = Msg("The input table contains no continuous spectra.")
        invalid_wavelength = Msg("The Wavelength meta contains missing values.")

    def __init__(self) -> None:
        super().__init__()
        self.data: Table | None = None
        self._build_controls()
        self._build_plot()

    def _build_controls(self) -> None:
        options = QGroupBox("Plot options", self.controlArea)
        layout = QVBoxLayout(options)

        layout.addWidget(QLabel("Processing stage", options))
        self.stage_combo = QComboBox(options)
        self.stage_combo.addItems(PROCESSING_STAGES)
        index = self.stage_combo.findText(self.selected_stage)
        self.stage_combo.setCurrentIndex(index if index >= 0 else 0)
        self.stage_combo.currentTextChanged.connect(self._stage_changed)
        layout.addWidget(self.stage_combo)

        layout.addWidget(QLabel("Colour scale", options))
        self.colour_scale_combo = QComboBox(options)
        self.colour_scale_combo.addItems(COLOUR_SCALES)
        index = self.colour_scale_combo.findText(self.colour_scale)
        self.colour_scale_combo.setCurrentIndex(index if index >= 0 else 0)
        self.colour_scale_combo.currentTextChanged.connect(
            self._colour_scale_changed
        )
        layout.addWidget(self.colour_scale_combo)

        layout.addWidget(QLabel("Line width", options))
        self.line_width_spin = QDoubleSpinBox(options)
        self.line_width_spin.setRange(0.1, 10.0)
        self.line_width_spin.setSingleStep(0.1)
        self.line_width_spin.setDecimals(1)
        self.line_width_spin.setSuffix(" px")
        self.line_width_spin.setValue(self.line_width)
        self.line_width_spin.valueChanged.connect(
            self._line_width_changed
        )
        layout.addWidget(self.line_width_spin)

        layout.addWidget(QLabel("Spectra", options))
        self.spectra_list = QListWidget(options)
        self.spectra_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.spectra_list.setMinimumWidth(280)
        self.spectra_list.setMinimumHeight(220)
        self.spectra_list.itemSelectionChanged.connect(self._replot)
        layout.addWidget(self.spectra_list)

        select_all = QPushButton("Select all", options)
        select_all.clicked.connect(self.spectra_list.selectAll)
        layout.addWidget(select_all)

        self.legend_checkbox = QCheckBox("Show legend", options)
        self.legend_checkbox.setChecked(self.show_legend)
        self.legend_checkbox.toggled.connect(self._legend_changed)
        layout.addWidget(self.legend_checkbox)

        self.reverse_checkbox = QCheckBox("Reverse wavelength axis", options)
        self.reverse_checkbox.setChecked(self.reverse_wavelength_axis)
        self.reverse_checkbox.toggled.connect(self._reverse_axis_changed)
        layout.addWidget(self.reverse_checkbox)

        self.controlArea.layout().addWidget(options)
        self.controlArea.layout().addStretch(1)

    def _build_plot(self) -> None:
        pg.setConfigOption("background", "w")
        pg.setConfigOption("foreground", "k")

        self.graphics_layout = pg.GraphicsLayoutWidget(self.mainArea)
        self.plot_item = self.graphics_layout.addPlot(row=0, col=0)
        self.plot_item.setLabel("bottom", "Wavelength")
        self.plot_item.setLabel("left", "Circular dichroism")
        self.plot_item.showGrid(x=True, y=True, alpha=0.2)

        # A separate layout column keeps the legend outside the plot frame.
        self.legend = pg.LegendItem(offset=(0, 0))
        self.graphics_layout.addItem(self.legend, row=0, col=1)
        self.graphics_layout.ci.layout.setColumnStretchFactor(0, 1)
        self.graphics_layout.ci.layout.setColumnStretchFactor(1, 0)
        self.mainArea.layout().addWidget(self.graphics_layout)

    @Inputs.data
    def set_data(self, data: Table | None) -> None:
        self.Error.clear()
        self.data = data
        self._populate_spectra()
        self._replot()

    def _stage_changed(self, stage: str) -> None:
        self.selected_stage = stage
        self._populate_spectra()
        self._replot()

    def _colour_scale_changed(self, scale_name: str) -> None:
        self.colour_scale = scale_name
        self._replot()

    def _legend_changed(self, checked: bool) -> None:
        self.show_legend = checked
        self._replot()

    def _reverse_axis_changed(self, checked: bool) -> None:
        self.reverse_wavelength_axis = checked
        self._replot()

    def _matching_variables(self) -> list[ContinuousVariable]:
        if self.data is None:
            return []
        variables = [
            variable
            for variable in self.data.domain.attributes
            if isinstance(variable, ContinuousVariable)
        ]
        stage = self.stage_combo.currentText()
        if stage == "All spectra":
            return variables
        suffix = f" | {stage}"
        return [variable for variable in variables if variable.name.endswith(suffix)]

    def _populate_spectra(self) -> None:
        selected = {
            item.data(Qt.UserRole) for item in self.spectra_list.selectedItems()
        }
        self.spectra_list.blockSignals(True)
        self.spectra_list.clear()
        for variable in self._matching_variables():
            item = QListWidgetItem(variable.name)
            item.setData(Qt.UserRole, variable.name)
            self.spectra_list.addItem(item)
            if not selected or variable.name in selected:
                item.setSelected(True)
        self.spectra_list.blockSignals(False)

    def _wavelength(self) -> np.ndarray | None:
        if self.data is None:
            return None
        try:
            variable = self.data.domain["Wavelength"]
        except KeyError:
            self.Error.missing_wavelength()
            return None
        if variable not in self.data.domain.metas or not isinstance(
            variable, ContinuousVariable
        ):
            self.Error.missing_wavelength()
            return None
        wavelength = np.asarray(self.data.get_column(variable), dtype=float)
        if np.isnan(wavelength).any():
            self.Error.invalid_wavelength()
            return None
        return wavelength

    def _line_width_changed(self, width: float) -> None:
        self.line_width = width
        self._replot()

    def _replot(self) -> None:
        self.plot_item.clear()
        self.legend.clear()
        self.plot_item.addLine(y=0, pen=pg.mkPen("#777777"))
        self.legend.setVisible(self.show_legend)
        self.Error.clear()

        if self.data is None:
            return
        wavelength = self._wavelength()
        if wavelength is None:
            return

        variable_by_name = {
            variable.name: variable for variable in self._matching_variables()
        }
        selected_names = [
            item.data(Qt.UserRole) for item in self.spectra_list.selectedItems()
        ]
        selected_variables = [
            variable_by_name[name]
            for name in selected_names
            if name in variable_by_name
        ]
        if not selected_variables:
            if not variable_by_name:
                self.Error.no_numeric_spectra()
            return

        colours = colours_from_scale(
            self.colour_scale_combo.currentText(), len(selected_variables)
        )
        for colour, variable in zip(colours, selected_variables):
            intensity = np.asarray(self.data.get_column(variable), dtype=float)
            valid = np.isfinite(wavelength) & np.isfinite(intensity)
            display_name = variable.name.split(" | ", maxsplit=1)[0]
            curve = self.plot_item.plot(
                wavelength[valid],
                intensity[valid],
                pen=pg.mkPen(colour, 
                             width=self.line_width_spin.value(),
                             ),
            )
            if self.show_legend:
                self.legend.addItem(curve, display_name)

        self.plot_item.invertX(self.reverse_wavelength_axis)
        self.plot_item.enableAutoRange()


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWCDSpectraPlot).run()