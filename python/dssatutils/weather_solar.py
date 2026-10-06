# Shared solar radiation helpers (extraterrestrial radiation and Bristow-Campbell estimation).
from __future__ import annotations

import math
import numpy as np
import pandas as pd


def extraterrestrial_radiation(lat_deg: float, doy: np.ndarray | list | float | int):
    """Daily extraterrestrial radiation Ra (MJ/m²/day) and sunset hour angle ws (rad).

    Implements FAO-56 Eq. 21 (Allen et al. 1998).
    Vectorised over day-of-year for a fixed latitude.
    Safely clips the sunset hour angle argument to [-1, 1] for polar day/night limits.
    """
    phi = math.radians(lat_deg)
    doy_arr = np.asarray(doy, dtype=float)
    dr = 1.0 + 0.033 * np.cos(2.0 * math.pi / 365.0 * doy_arr)          # inverse relative distance
    decl = 0.409 * np.sin(2.0 * math.pi / 365.0 * doy_arr - 1.39)       # solar declination
    arg = np.clip(-np.tan(phi) * np.tan(decl), -1.0, 1.0)
    ws = np.arccos(arg)                                                 # sunset hour angle
    Gsc = 0.0820                                                        # solar constant (MJ/m²/min)
    Ra = (24.0 * 60.0 / math.pi) * Gsc * dr * (
        ws * math.sin(phi) * np.sin(decl)
        + math.cos(phi) * np.cos(decl) * np.sin(ws)
    )
    return np.maximum(Ra, 0.0), ws


def estimate_srad_bristow_campbell(
    dates,
    tmax,
    tmin,
    lat_deg: float,
    a: float = 0.75,
    c: float = 2.4,
) -> np.ndarray:
    """Estimate daily solar radiation (MJ/m²/day) from temperature range.

    Method: Bristow, K. L., & Campbell, G. S. (1984). On the relationship between
    incoming solar radiation and daily maximum and minimum temperature.
    Agricultural and Forest Meteorology, 31(2), 159-166.

    Parameters
    ----------
    dates : sequence of dates
        Calendar dates (DatetimeIndex, Series of dates/strings/YYYYDOY, etc.).
    tmax : array-like
        Daily maximum air temperature (°C).
    tmin : array-like
        Daily minimum air temperature (°C).
    lat_deg : float
        Latitude in decimal degrees (-90 to +90).
    a : float, default 0.75
        Maximum atmospheric transmittance (clear-sky coefficient).
    c : float, default 2.4
        Empirical shape exponent.

    Returns
    -------
    np.ndarray
        Estimated daily solar radiation (MJ/m²/day), clipped to [0, 0.80 * Ra].
        NaN where Tmax or Tmin are missing/non-finite.
    """
    tmax_arr = np.asarray(tmax, dtype=float)
    tmin_arr = np.asarray(tmin, dtype=float)
    n = len(tmax_arr)
    if n == 0:
        return np.array([], dtype=float)

    # Parse dates to get DOY, calendar month, and check day-step continuity
    if isinstance(dates, pd.DatetimeIndex):
        dt = dates
    elif isinstance(dates, pd.Series) and pd.api.types.is_datetime64_any_dtype(dates):
        dt = pd.DatetimeIndex(dates)
    else:
        s = pd.Series(dates).astype(str)
        if len(s) > 0 and s.iloc[0].isdigit() and len(s.iloc[0]) == 7:
            dt = pd.DatetimeIndex(pd.to_datetime(s, format="%Y%j"))
        else:
            dt = pd.DatetimeIndex(pd.to_datetime(s))

    doy = dt.dayofyear.values
    months = dt.month.values

    # Diurnal temperature range: Delta_T(t) = Tmax(t) - [Tmin(t) + Tmin(t+1)] / 2
    # For the terminal day or across date gaps (non-consecutive days), Delta_T(t) = Tmax(t) - Tmin(t).
    tmin_next = np.empty(n, dtype=float)
    tmin_next[:-1] = tmin_arr[1:]
    tmin_next[-1] = tmin_arr[-1]

    # Flag whether day t+1 is exactly 1 calendar day after day t
    if n > 1:
        consecutive = np.empty(n, dtype=bool)
        consecutive[:-1] = (dt[1:].values - dt[:-1].values) == np.timedelta64(1, "D")
        consecutive[-1] = False
    else:
        consecutive = np.array([False], dtype=bool)

    use_tmin_next = consecutive & np.isfinite(tmin_next)
    tmin_effective = np.where(use_tmin_next, (tmin_arr + tmin_next) / 2.0, tmin_arr)
    delta_t = tmax_arr - tmin_effective
    delta_t = np.where(np.isfinite(delta_t), np.maximum(delta_t, 0.0), np.nan)

    # Monthly mean diurnal range Delta_T_bar_m across valid days
    valid_dt = delta_t[np.isfinite(delta_t) & (delta_t > 0)]
    fallback_dt_bar = float(np.mean(valid_dt)) if len(valid_dt) > 0 else 10.0
    if not np.isfinite(fallback_dt_bar) or fallback_dt_bar <= 0:
        fallback_dt_bar = 10.0

    dt_bar_m = {}
    for m in range(1, 13):
        m_mask = (months == m) & np.isfinite(delta_t) & (delta_t > 0)
        if np.any(m_mask):
            dt_bar_m[m] = float(np.mean(delta_t[m_mask]))
        else:
            dt_bar_m[m] = fallback_dt_bar

    dt_bar_series = np.array([dt_bar_m[m] for m in months], dtype=float)
    b = 0.036 * np.exp(-0.154 * dt_bar_series)

    # Atmospheric transmittance Tt = A * [1 - exp(-B * Delta_T^C)]
    tt = a * (1.0 - np.exp(-b * (delta_t ** c)))

    # Extraterrestrial radiation Ra
    ra, _ = extraterrestrial_radiation(lat_deg, doy)

    # Solar radiation Rs = Tt * Ra clipped to [0, 0.80 * Ra]
    rs = tt * ra
    rs = np.where(ra <= 0.0, 0.0, rs)
    rs = np.clip(rs, 0.0, 0.80 * ra)

    # Missing inputs produce missing output
    invalid = ~np.isfinite(tmax_arr) | ~np.isfinite(tmin_arr)
    rs = np.where(invalid, np.nan, rs)
    return rs
