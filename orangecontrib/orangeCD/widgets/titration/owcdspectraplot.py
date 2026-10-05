#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Orange widget for plotting circular-dichroism spectra."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from AnyQt.QtGui import QColor
from AnyQt.QtWidgets import QAbstractItemView
from Orange.data import Table
from AnyQt.QtCore import Qt
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import Input, Msg, OWWidget

from .utils import (
    InvalidWavelength,
    SpectraError,
    matching_spectra,
    WAVELENGTH_UNIT_KEY,
    shared_unit,
    spectrum_names,
    spectrum_units,
    table_unit,
    unit_symbol,
    wavelengths as spectra_wavelengths,
)

COLOUR_SCALES = {
    "Viridis":      ("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"),
    "Plasma":       ("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921"),
    "Inferno":      ("#000004", "#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4"),
    "Cividis":      ("#00224e", "#31446b", "#666970", "#a38f63", "#e6c75a", "#fee838"),
    "Blue to red":  ("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b"),
    "Greyscale":    ("#111111", "#555555", "#999999", "#dddddd"),
    "Classic":      ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
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
        first = anchors[lower]
        second = anchors[upper]
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
    description = "Plot spectra produced by the CD Data Loader widget or by spectral preprocessing."
    icon = "icons/CDSpectraPlot.svg"
    priority = 30
    want_main_area = True
    resizing_enabled = True

    selected_stage = Setting("All spectra")
    colour_scale = Setting("Viridis")
    line_width = Setting(1.5)
    show_legend = Setting(True)
    reverse_wavelength_axis = Setting(False)
    selected_spectra = Setting([])
    processing_stages = Setting(["All spectra"])

    # This is the list model displayed by gui.listBox. It is derived from the
    # input domain, so it is not persisted as a Setting.
    spectra_names: list[str] = []

    class Inputs:
        data = Input("CD Data", Table)

    class Error(OWWidget.Error):
        missing_wavelength = Msg("The input table does not contain a 'Spectrum' string meta naming the rows.")
        no_numeric_spectra = Msg("The input table contains no continuous spectra.")
        invalid_wavelength = Msg("The input attribute names are not all valid wavelengths.")

    def __init__(self) -> None:
        super().__init__()

        self.colour_scale = self._normalise_combo_setting(
            self.colour_scale,
            tuple(COLOUR_SCALES),
            default="Viridis",
        )

        self.data: Table | None = None
        # Maps list-box position -> row index in self.data
        self._row_map: list[int] = []
        self._build_controls()
        self._build_plot()
        
    def _build_controls(self) -> None:
        options = gui.widgetBox(self.controlArea, "Plot options")

        self.stage_combo = gui.comboBox(
            options,
            self,
            "selected_stage",
            label="Processing stage",
            items=self.processing_stages,
            orientation=Qt.Vertical,
            sendSelectedValue=True,
            callback=self._stage_changed,
        )

        gui.comboBox(
            options,
            self,
            "colour_scale",
            label="Colour scale",
            items=tuple(COLOUR_SCALES),
            orientation=Qt.Vertical,
            sendSelectedValue=True,
            callback=self._replot,
        )

        gui.doubleSpin(
            options,
            self,
            "line_width",
            0.1,
            10.0,
            step=0.1,
            label="Line width",
            decimals=1,
            suffix=" px",
            orientation=Qt.Vertical,
            callback=self._replot,
        )

        gui.widgetLabel(options, "Spectra")
        self.spectra_list = gui.listBox(
            options,
            self,
            "selected_spectra",
            "spectra_names",
            selectionMode=QAbstractItemView.ExtendedSelection,
            callback=self._replot,
        )
        self.spectra_list.setMinimumWidth(280)
        self.spectra_list.setMinimumHeight(220)

        gui.button(
            options,
            self,
            "Select all",
            callback=self._select_all_spectra,
        )

        gui.checkBox(
            options,
            self,
            "show_legend",
            "Show legend",
            callback=self._replot,
        )
        gui.checkBox(
            options,
            self,
            "reverse_wavelength_axis",
            "Reverse wavelength axis",
            callback=self._replot,
        )
        gui.rubber(self.controlArea)

    @staticmethod
    def _normalise_combo_setting(
        value: str | int,
        choices: tuple[str, ...],
        *,
        default: str,
    ) -> str:
        """Convert settings saved by an index-based combo box to text."""

        if isinstance(value, int):
            if 0 <= value < len(choices):
                return choices[value]
            return default

        if value in choices:
            return value

        return default
        
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

        self._set_processing_stages()
        self._populate_spectra()
        self._replot()

    def _set_processing_stages(self) -> None:
        current = self.selected_stage

        if self.data is None:
            stages = ["All spectra"]
        else:
            try:
                names = spectrum_names(self.data)
            except SpectraError:
                names = []
            stages = sorted({name.split("|")[-1].strip() for name in names})
            stages.insert(0, "All spectra")
        self.processing_stages = stages

        self.stage_combo.clear()
        self.stage_combo.addItems(stages)

        if current in stages:
            self.selected_stage = current
            self.stage_combo.setCurrentText(current)
        else:
            self.selected_stage = stages[0]
            self.stage_combo.setCurrentIndex(0)
        
    def _stage_changed(self) -> None:
        self._populate_spectra()
        self._replot()

    def _matching_variables(self) -> list[tuple[int, str]]:
        if self.data is None:
            return []
        stage = None if self.selected_stage == "All spectra" else self.selected_stage
        try:
            return matching_spectra(self.data, stage)
        except SpectraError:
            return []

    def _populate_spectra(self) -> None:
        matched_vars = self._matching_variables()

        self._row_map = [i[0] for i in matched_vars]
        self.spectra_names = [i[1] for i in matched_vars]
        self.selected_spectra = list(range(len(matched_vars)))

    def _select_all_spectra(self) -> None:
        self.selected_spectra = list(range(len(self.spectra_names)))
        self._replot()

    def _update_axis_labels(self) -> None:
        wavelength_unit = None
        signal_unit = None

        if self.data is not None:
            wavelength_unit = table_unit(self.data, WAVELENGTH_UNIT_KEY)

            # Only label the y axis when all plotted spectra share one unit.
            row_units = spectrum_units(self.data)
            plotted = [
                row_units[self._row_map[index]]
                for index in self.selected_spectra
                if 0 <= index < len(self._row_map)
            ]
            if plotted:
                signal_unit = shared_unit(plotted)

        self.plot_item.setLabel(
            "bottom", "Wavelength", units=unit_symbol(wavelength_unit)
        )
        self.plot_item.setLabel(
            "left", "Circular dichroism", units=unit_symbol(signal_unit)
        )

    def _wavelength(self) -> np.ndarray | None:
        if self.data is None:
            return None

        try:
            return spectra_wavelengths(self.data)
        except InvalidWavelength:
            self.Error.invalid_wavelength()
        except SpectraError:
            self.Error.missing_wavelength()
        return None

    def _replot(self) -> None:
        self.plot_item.clear()
        self.legend.clear()
        self.plot_item.addLine(y=0, pen=pg.mkPen("#777777"))
        self.legend.setVisible(self.show_legend)
        self.Error.clear()
        self._update_axis_labels()

        if self.data is None:
            return

        wavelength = self._wavelength()
        if wavelength is None:
            return

        selected = [
            index for index in self.selected_spectra
            if 0 <= index < len(self._row_map)
        ]
        if not selected:
            if not self.spectra_names:
                self.Error.no_numeric_spectra()
            return

        colours = colours_from_scale(
            self.colour_scale,
            len(selected),
        )
        for colour, index in zip(colours, selected):
            var_name = self.spectra_names[index]
            intensity = self.data.X[self._row_map[index]]
            valid = np.isfinite(wavelength) & np.isfinite(intensity)
            curve = self.plot_item.plot(
                wavelength[valid],
                intensity[valid],
                pen=pg.mkPen(colour, width=self.line_width),
            )
            if self.show_legend:
                display_name = var_name.split(" | ", maxsplit=1)[0]
                self.legend.addItem(curve, display_name)

        self.plot_item.invertX(self.reverse_wavelength_axis)
        self.plot_item.enableAutoRange()


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWCDSpectraPlot).run()
