#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Aug  6 17:41:17 2026

@author: ubx84221
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Sequence

import pandas as pd


class TitrationMode(str, Enum):
    FIXED = "fixed"
    INCREASING = "increasing"


@dataclass
class TitrationPoint:
    ratio: float
    stock_b: int
    predicted_volume: float | None = None


@dataclass
class IncreasingRow:
    ratio: float
    stock_b: int

    predicted_volume: float
    volume_added_this_step: float
    total_stock_b_volume: float
    total_cell_volume: float
    dilution_factor: float
    normalised_molar_ratio: float


@dataclass
class FixedRow:
    ratio: float
    stock_b: int

    predicted_volume: float
    volume_stock_b: float
    baseline_volume: float
    concentration_b: float
    cell_volume: float
    normalised_molar_ratio: float


@dataclass
class TitrationResult:
    mode: TitrationMode
    volume_solution_a: float

    # Increasing-volume mode only
    volume_buffer: float | None = None
    max_volume_allowed: float | None = None
    max_volume_added: float | None = None
    within_limit: bool | None = None

    rows: list[IncreasingRow] | list[FixedRow] = field(
        default_factory=list
    )

    def result_to_dataframe(self) -> pd.DataFrame:
        """
        Convert the calculated titration rows into a pandas DataFrame.

        Calculation-level metadata is added as additional columns.
        """

        df = pd.DataFrame(asdict(row) for row in self.rows)

        df["mode"] = self.mode.value
        df["volume_solution_a"] = self.volume_solution_a

        if self.mode == TitrationMode.INCREASING:
            df["volume_buffer"] = self.volume_buffer
            df["max_volume_allowed"] = self.max_volume_allowed
            df["max_volume_added"] = self.max_volume_added
            df["within_limit"] = self.within_limit

        return df


class CDTitrationCalculator:
    """
    Reimplementation of the frmTitration calculation logic.

    Concentrations must all use the same concentration unit.
    Volumes are returned in the same volume unit as
    starting_cell_volume.
    """

    def __init__(
        self,
        starting_cell_volume: float,
        stock_con_a: float,
        working_con_a: float,
        stock_b_concentrations: Sequence[float],
        stock_b_molar_equiv: float,
    ) -> None:
        if starting_cell_volume <= 0:
            raise ValueError(
                "Starting cell volume must be greater than zero"
            )

        if stock_con_a <= 0:
            raise ValueError(
                "Stock concentration of solution A must be "
                "greater than zero"
            )

        if working_con_a <= 0:
            raise ValueError(
                "Working concentration of solution A must be "
                "greater than zero"
            )

        if not stock_b_concentrations:
            raise ValueError(
                "At least one stock-B concentration is required"
            )

        if any(
            concentration <= 0
            for concentration in stock_b_concentrations
        ):
            raise ValueError(
                "All stock-B concentrations must be greater "
                "than zero"
            )

        self.starting_cell_volume = float(
            starting_cell_volume
        )
        self.stock_con_a = float(stock_con_a)
        self.working_con_a = float(working_con_a)
        self.stock_b_concs = [
            float(concentration)
            for concentration in stock_b_concentrations
        ]
        self.stock_b_molar_equiv = float(stock_b_molar_equiv)

    @property
    def volume_solution_a(self) -> float:
        """
        Constant volume of stock solution A required.
        """

        return round(
            (
                self.working_con_a
                / self.stock_con_a
            )
            * self.starting_cell_volume,
            1,
        )

    def _required_stock_b_volume(
        self,
        required_ratio: float,
        stock_con_b: float,
    ) -> float:
        """
        Calculate the unrounded stock-B volume required for a
        ratio, or for a ratio increment in increasing mode.
        """

        if required_ratio < 0:
            raise ValueError(
                "The required ratio must not be negative"
            )

        return (
            required_ratio
            * self.working_con_a
            * self.starting_cell_volume
            / stock_con_b
        )

    @staticmethod
    def _distance_from_target_range(
        volume: float,
        target_min: float,
        target_max: float,
    ) -> float:
        """
        Return the distance from a volume to a closed target range.

        A volume inside the target range has a distance of zero.
        """

        if volume < target_min:
            return target_min - volume

        if volume > target_max:
            return volume - target_max

        return 0.0

    def choose_stock(
        self,
        required_ratio: float,
        *,
        target_min: float = 2.0,
        target_max: float = 20.0,
    ) -> tuple[int, float]:
        """
        Automatically choose the most appropriate stock-B solution.

        The least concentrated stock whose required volume falls
        inside the target pipetting range is preferred. Because
        stock_b_concs retains its supplied order, it should normally
        be ordered from the least concentrated stock to the most
        concentrated stock.

        If no stock produces a volume inside the target range, the
        stock producing the volume closest to that range is selected.

        Parameters
        ----------
        required_ratio
            In fixed-volume mode this is the full B:A molar ratio.

            In increasing-volume mode this is the increment in B:A
            molar ratio since the preceding titration point.

        target_min
            Preferred minimum pipetting volume.

        target_max
            Preferred maximum pipetting volume.

        Returns
        -------
        tuple[int, float]
            Stock number and unrounded predicted pipetting volume.
        """

        if required_ratio <= 0:
            raise ValueError(
                "Required ratio must be greater than zero"
            )

        if target_min < 0:
            raise ValueError(
                "Minimum target volume must not be negative"
            )

        if target_max <= target_min:
            raise ValueError(
                "Maximum target volume must be greater than "
                "minimum target volume"
            )

        candidates: list[tuple[int, float]] = []

        for stock_number, stock_con_b in enumerate(
            self.stock_b_concs,
            start=1,
        ):
            predicted_volume = (
                self._required_stock_b_volume(
                    required_ratio=required_ratio,
                    stock_con_b=stock_con_b,
                )
            )

            candidates.append(
                (stock_number, predicted_volume)
            )

        # Prefer the first, normally least concentrated, stock
        # producing a volume within the desired pipetting range.
        for stock_number, predicted_volume in candidates:
            if target_min <= predicted_volume <= target_max:
                return stock_number, predicted_volume

        # No candidate is inside the target range. Choose the one
        # whose calculated volume is closest to the range.
        return min(
            candidates,
            key=lambda candidate: (
                self._distance_from_target_range(
                    volume=candidate[1],
                    target_min=target_min,
                    target_max=target_max,
                ),
                candidate[0],
            ),
        )

    def create_points(
        self,
        ratios: Sequence[float],
        *,
        mode: TitrationMode,
        target_min: float = 2.0,
        target_max: float = 20.0,
    ) -> list[TitrationPoint]:
        """
        Create titration points using automatic stock selection.

        For fixed-volume titrations, each stock is selected using the
        complete ratio at that point.

        For increasing-volume titrations, each stock is selected
        using the ratio increment since the preceding point, because
        that increment determines the volume added at that step.
        """

        numeric_ratios = [
            float(ratio)
            for ratio in ratios
        ]

        if not numeric_ratios:
            raise ValueError(
                "At least one titration ratio is required"
            )

        if any(ratio <= 0 for ratio in numeric_ratios):
            raise ValueError(
                "All titration ratios must be greater than zero"
            )

        if mode == TitrationMode.INCREASING:
            if any(
                current <= previous
                for previous, current in zip(
                    numeric_ratios,
                    numeric_ratios[1:],
                )
            ):
                raise ValueError(
                    "Increasing-mode titration ratios must be "
                    "strictly increasing"
                )

        points: list[TitrationPoint] = []
        previous_ratio = 0.0

        for ratio in numeric_ratios:
            if mode == TitrationMode.FIXED:
                required_ratio = ratio
            elif mode == TitrationMode.INCREASING:
                required_ratio = ratio - previous_ratio
            else:
                raise ValueError(
                    f"Unknown titration mode: {mode}"
                )

            stock_b, predicted_volume = self.choose_stock(
                required_ratio=required_ratio,
                target_min=target_min,
                target_max=target_max,
            )

            points.append(
                TitrationPoint(
                    ratio=ratio,
                    stock_b=stock_b,
                    predicted_volume=predicted_volume,
                )
            )

            previous_ratio = ratio

        return points

    def _validate_points(
        self,
        points: Sequence[TitrationPoint],
    ) -> None:
        """
        Validate manually or automatically generated points.
        """

        if not points:
            raise ValueError(
                "At least one titration point is required"
            )

        number_of_stocks = len(self.stock_b_concs)

        for point in points:
            if point.ratio <= 0:
                raise ValueError(
                    "All titration ratios must be greater than zero"
                )

            if not 1 <= point.stock_b <= number_of_stocks:
                raise ValueError(
                    f"Stock-B number {point.stock_b} is invalid. "
                    f"Expected a value from 1 to "
                    f"{number_of_stocks}."
                )

    def calculate(
        self,
        mode: TitrationMode,
        points: Sequence[TitrationPoint],
    ) -> TitrationResult:
        """
        Calculate the requested titration table.
        """

        self._validate_points(points)

        if mode == TitrationMode.FIXED:
            return self._calculate_fixed(points)

        if mode == TitrationMode.INCREASING:
            return self._calculate_increasing(points)

        raise ValueError(f"Unknown mode: {mode}")

    def _calculate_fixed(
        self,
        points: Sequence[TitrationPoint],
    ) -> TitrationResult:
        """
        Calculate a fixed-final-volume titration.
        """

        volume_a = self.volume_solution_a
        rows: list[FixedRow] = []

        for point in points:
            stock_b_con = self.stock_b_concs[
                point.stock_b - 1
            ]

            calculated_volume = (
                self._required_stock_b_volume(
                    required_ratio=point.ratio,
                    stock_con_b=stock_b_con,
                )
            )

            predicted_volume = (
                calculated_volume
                if point.predicted_volume is None
                else point.predicted_volume
            )

            volume_stock_b = round(
                calculated_volume,
                1,
            )

            baseline_volume = round(
                self.starting_cell_volume
                - volume_a
                - volume_stock_b,
                1,
            )

            if baseline_volume < 0:
                raise ValueError(
                    f"The calculated baseline volume is negative "
                    f"at ratio {point.ratio}. The combined volumes "
                    f"of solutions A and B exceed the fixed cell "
                    f"volume."
                )

            concentration_b = round(
                self.working_con_a * point.ratio,
                1,
            )

            rows.append(
                FixedRow(
                    ratio=point.ratio,
                    stock_b=point.stock_b,
                    predicted_volume=round(
                        predicted_volume,
                        3,
                    ),
                    volume_stock_b=volume_stock_b,
                    baseline_volume=baseline_volume,
                    concentration_b=concentration_b,
                    cell_volume=self.starting_cell_volume,
                    normalised_molar_ratio=round(point.ratio/self.stock_b_molar_equiv,3)
                )
            )

        return TitrationResult(
            mode=TitrationMode.FIXED,
            volume_solution_a=volume_a,
            rows=rows,
        )

    def _calculate_increasing(
        self,
        points: Sequence[TitrationPoint],
    ) -> TitrationResult:
        """
        Calculate an increasing-volume titration.
        """

        volume_a = self.volume_solution_a

        volume_buffer = round(
            self.starting_cell_volume - volume_a,
            1,
        )

        rows: list[IncreasingRow] = []

        previous_ratio = 0.0
        total_stock_b = 0.0

        for point in points:
            if point.ratio <= previous_ratio:
                raise ValueError(
                    "Titration points must be strictly increasing"
                )

            stock_b_con = self.stock_b_concs[
                point.stock_b - 1
            ]

            ratio_increment = (
                point.ratio - previous_ratio
            )

            calculated_step_volume = (
                self._required_stock_b_volume(
                    required_ratio=ratio_increment,
                    stock_con_b=stock_b_con,
                )
            )

            predicted_volume = (
                calculated_step_volume
                if point.predicted_volume is None
                else point.predicted_volume
            )

            step_volume = round(
                calculated_step_volume,
                1,
            )

            total_stock_b = round(
                total_stock_b + step_volume,
                1,
            )

            total_cell_volume = round(
                self.starting_cell_volume
                + total_stock_b,
                1,
            )

            dilution_factor = round(
                total_cell_volume
                / self.starting_cell_volume,
                3,
            )

            rows.append(
                IncreasingRow(
                    ratio=point.ratio,
                    stock_b=point.stock_b,
                    predicted_volume=round(
                        predicted_volume,
                        3,
                    ),
                    volume_added_this_step=step_volume,
                    total_stock_b_volume=total_stock_b,
                    total_cell_volume=total_cell_volume,
                    dilution_factor=dilution_factor,
                    normalised_molar_ratio=round(point.ratio/self.stock_b_molar_equiv,3)
                )
            )

            previous_ratio = point.ratio

        max_volume_allowed = round(
            self.starting_cell_volume * 0.15,
            1,
        )

        return TitrationResult(
            mode=TitrationMode.INCREASING,
            volume_solution_a=volume_a,
            volume_buffer=volume_buffer,
            max_volume_allowed=max_volume_allowed,
            max_volume_added=total_stock_b,
            within_limit=(
                total_stock_b <= max_volume_allowed
            ),
            rows=rows,
        )


if __name__ == "__main__":
    calculator = CDTitrationCalculator(
        starting_cell_volume=500.0,
        stock_con_a=468.0,
        working_con_a=19.659,
        stock_b_concentrations=[
            2000.0,
            4000.0,
            8000.0,
        ],
        stock_b_molar_equiv=2.75
    )

    ratios = [
        0.14,
        0.28,
        0.42,
        0.56,
        0.70,
        0.84,
        0.98,
        1.12,
        1.26,
        1.40,
        1.68,
        1.96,
        2.24,
        2.52,
        2.80,
        10,
        30
    ]

    target_min_volume = 2.0
    target_max_volume = 20.0

    fixed_points = calculator.create_points(
        ratios=ratios,
        mode=TitrationMode.FIXED,
        target_min=target_min_volume,
        target_max=target_max_volume,
    )

    fixed_result = calculator.calculate(
        mode=TitrationMode.FIXED,
        points=fixed_points,
    )

    print("Fixed-volume mode")
    print(
        "Volume of solution A:",
        fixed_result.volume_solution_a,
    )
    print(fixed_result.result_to_dataframe().to_string(index=False))

    increasing_points = calculator.create_points(
        ratios=ratios,
        mode=TitrationMode.INCREASING,
        target_min=target_min_volume,
        target_max=target_max_volume,
    )

    increasing_result = calculator.calculate(
        mode=TitrationMode.INCREASING,
        points=increasing_points,
    )

    print("\nIncreasing-volume mode")
    print(
        "Volume of solution A:",
        increasing_result.volume_solution_a,
    )
    print(
        "Initial buffer volume:",
        increasing_result.volume_buffer,
    )
    print(
        increasing_result
        .result_to_dataframe()
        .to_string(index=False)
    )