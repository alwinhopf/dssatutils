"""Helpers for replacing coarse weather-source rainfall with gridded rain."""

from __future__ import annotations

from typing import Mapping

import pandas as pd


def merge_rainfall_into_weather(
    weather: pd.DataFrame,
    rainfall: Mapping[str, float] | pd.Series | None,
    date_col: str = "DATE",
    rain_col: str = "RAIN",
) -> int:
    """Replace *weather[rain_col]* where daily rainfall has matching DATE keys.

    Parameters
    ----------
    weather:
        Daily weather frame with DSSAT-style ``YYYYDOY`` date codes.
    rainfall:
        Mapping or Series keyed by ``YYYYDOY`` with rainfall in mm/day.
    date_col, rain_col:
        Column names in *weather*.

    Returns
    -------
    int
        Number of daily rows replaced. The input DataFrame is modified in place.
    """
    if rainfall is None:
        return 0
    series = pd.Series(rainfall, dtype="float64")
    if series.empty:
        return 0
    series.index = series.index.astype(str)
    mapped = weather[date_col].astype(str).map(series)
    use = mapped.notna()
    if use.any():
        weather[rain_col] = pd.to_numeric(weather[rain_col], errors="coerce").astype("float64")
        weather.loc[use, rain_col] = mapped[use].values
    return int(use.sum())
