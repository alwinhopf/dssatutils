# File: weather_dwd.py
# Python port of weather_dwd.R
#
# Weather source: DWD (Deutscher Wetterdienst) Open Data — daily climate station
# observations for Germany -> DSSAT .WTH.
#
# WHY: for German / Central-European sites (e.g. the winter-wheat phenology work)
# the national met service's quality-controlled station network is more accurate
# than a global reanalysis. DWD's daily "kl" product gives max/min/mean air
# temperature, precipitation, sunshine duration, wind, humidity and vapour
# pressure directly from gauges.
#
# ACCESS (fully open, no key): https://opendata.dwd.de/ (CDC).
#   Stations:  .../climate/daily/kl/historical/KL_Tageswerte_Beschreibung_Stationen.txt
#   Per-station daily:  historical/tageswerte_KL_<id>_<from>_<to>_hist.zip  (+ recent akt zip)
#
# DSSAT needs daily SOLAR RADIATION, which the kl product does NOT measure
# directly; it is estimated from sunshine duration (SDK, hours) with the
# Angstrom-Prescott relation (Rs = (a_s + b_s n/N) Ra, FAO-56), which is the
# standard fallback when measured Rs is unavailable. Where SDK is missing, SRAD
# is written as -99 (DSSAT-valid missing).
#
# Coverage: Germany. The network fetch (_dwd_stations / _fetch_station) is
# isolated from the .WTH formatting (_write_wth) so the latter is unit-testable
# with synthetic data, exactly like the AgERA5 module.

import io
import math
import os
import zipfile
from datetime import date

import numpy as np
import pandas as pd
import requests

_DWD_KL = ("https://opendata.dwd.de/climate_environment/CDC/"
           "observations_germany/climate/daily/kl/")
_STATION_DESC = "KL_Tageswerte_Beschreibung_Stationen.txt"

# DWD product column -> DSSAT variable. -999 is the DWD missing sentinel.
#   TXK air-temp max (°C), TNK min (°C), RSK precip (mm), SDK sunshine (h),
#   FM mean wind (m/s, 10 m), UPM relative humidity (%), VPM vapour pressure (hPa).


# ---------------------------------------------------------------------------
# Climatology helpers (shared convention with the other weather modules)
# ---------------------------------------------------------------------------

def _calc_tav(df: pd.DataFrame) -> float:
    return float(((df["TMAX"] + df["TMIN"]) / 2.0).mean())


def _calc_amp(df: pd.DataFrame) -> float:
    d = df.copy()
    d["TAVG"] = (d["TMAX"] + d["TMIN"]) / 2.0
    monthly = d.groupby(["YEAR", "MM"])["TAVG"].mean().reset_index()
    annual = monthly.groupby("YEAR")["TAVG"].agg(lambda x: x.max() - x.min())
    return float(annual.mean())


# ---------------------------------------------------------------------------
# Solar radiation from sunshine duration (Angstrom-Prescott / FAO-56)
# ---------------------------------------------------------------------------

def _extraterrestrial_radiation(lat_deg, doy):
    """Daily extraterrestrial radiation Ra (MJ/m²/day). FAO-56 eq. 21.

    Vectorised over arrays of day-of-year for a fixed latitude.
    """
    phi = math.radians(lat_deg)
    doy = np.asarray(doy, dtype=float)
    dr = 1.0 + 0.033 * np.cos(2.0 * math.pi / 365.0 * doy)            # inverse rel. distance
    decl = 0.409 * np.sin(2.0 * math.pi / 365.0 * doy - 1.39)         # solar declination
    arg = np.clip(-np.tan(phi) * np.tan(decl), -1.0, 1.0)
    ws = np.arccos(arg)                                               # sunset hour angle
    Gsc = 0.0820                                                      # MJ/m²/min
    Ra = (24.0 * 60.0 / math.pi) * Gsc * dr * (
        ws * math.sin(phi) * np.sin(decl)
        + math.cos(phi) * np.cos(decl) * np.sin(ws))
    return np.maximum(Ra, 0.0), ws


def _srad_from_sunshine(lat_deg, doy, sunshine_h, a_s=0.25, b_s=0.50):
    """Estimate daily solar radiation (MJ/m²/day) from sunshine hours.

    Returns NaN where sunshine is missing/NaN. Clipped to [0, 0.8*Ra] (clear-sky
    transmissivity ceiling) so cloud-free days don't overshoot.
    """
    Ra, ws = _extraterrestrial_radiation(lat_deg, doy)
    N = 24.0 / math.pi * ws                                          # max daylight hours
    n = np.asarray(sunshine_h, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(N > 0, np.clip(n / N, 0.0, 1.0), 0.0)
        Rs = (a_s + b_s * frac) * Ra
    Rs = np.where(np.isnan(n), np.nan, np.clip(Rs, 0.0, 0.8 * Ra))
    return Rs


def _tdew_from_vapour_pressure(vpm_hpa):
    """Dew-point (°C) from vapour pressure (hPa) via the inverse Magnus formula."""
    e = np.asarray(vpm_hpa, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        ln = np.log(np.where(e > 0, e / 6.1094, np.nan))
        td = (243.04 * ln) / (17.625 - ln)
    return td


# ---------------------------------------------------------------------------
# Station metadata + per-station daily data (network; cached on disk)
# ---------------------------------------------------------------------------

def _dwd_stations(cache_dir: str) -> pd.DataFrame:
    """Download + parse the DWD kl station description (id, lat, lon, elev, dates)."""
    os.makedirs(cache_dir, exist_ok=True)
    local = os.path.join(cache_dir, _STATION_DESC)
    if not (os.path.exists(local) and os.path.getsize(local) > 0):
        r = requests.get(_DWD_KL + "historical/" + _STATION_DESC, timeout=120)
        r.raise_for_status()
        with open(local, "wb") as fh:
            fh.write(r.content)
    with open(local, "r", encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    rows = []
    for ln in lines[2:]:                       # skip header + dashed rule
        if not ln.strip():
            continue
        # The first 6 fields are whitespace-delimited numbers; only the trailing
        # station name + federal state contain spaces, so split off 6 tokens.
        parts = ln.split(None, 6)
        if len(parts) < 6:
            continue
        try:
            rows.append({
                "station_id": parts[0].zfill(5),
                "von": int(parts[1]), "bis": int(parts[2]),
                "elev": float(parts[3]) if parts[3] else np.nan,
                "lat": float(parts[4]), "lon": float(parts[5]),
            })
        except Exception:  # noqa: BLE001
            continue
    if not rows:
        return pd.DataFrame(columns=["station_id", "von", "bis", "elev", "lat", "lon"])
    return pd.DataFrame(rows).dropna(subset=["lat", "lon"])


def _historical_index(cache_dir: str) -> dict:
    """Map station_id -> historical zip filename by listing the historical dir."""
    local = os.path.join(cache_dir, "_hist_index.txt")
    if os.path.exists(local) and os.path.getsize(local) > 0:
        names = open(local, encoding="utf-8").read().splitlines()
    else:
        import re
        r = requests.get(_DWD_KL + "historical/", timeout=120)
        r.raise_for_status()
        names = re.findall(r'href="(tageswerte_KL_[^"]+_hist\.zip)"', r.text)
        with open(local, "w", encoding="utf-8") as fh:
            fh.write("\n".join(names))
    idx = {}
    for n in names:
        parts = n.split("_")
        if len(parts) >= 3:
            idx[parts[2].zfill(5)] = n
    return idx


def _parse_product_zip(content: bytes) -> pd.DataFrame:
    """Parse a DWD kl product zip into a daily DataFrame keyed by date."""
    zf = zipfile.ZipFile(io.BytesIO(content))
    prod = next((n for n in zf.namelist() if n.startswith("produkt")), None)
    if prod is None:
        return pd.DataFrame()
    raw = zf.read(prod).decode("latin-1")
    df = pd.read_csv(io.StringIO(raw), sep=";", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    if "MESS_DATUM" not in df.columns:
        return pd.DataFrame()
    df["DATE"] = pd.to_datetime(df["MESS_DATUM"].astype(str), format="%Y%m%d", errors="coerce")
    df = df.dropna(subset=["DATE"])
    df = df.replace(-999, np.nan).replace(-999.0, np.nan)
    return df


def _fetch_station(station_id: str, hist_name, cache_dir: str) -> pd.DataFrame:
    """Download (historical + recent) and parse one station's daily series."""
    frames = []
    if hist_name:
        local = os.path.join(cache_dir, hist_name)
        if not (os.path.exists(local) and os.path.getsize(local) > 0):
            try:
                rr = requests.get(_DWD_KL + "historical/" + hist_name, timeout=180)
                if rr.ok:
                    open(local, "wb").write(rr.content)
            except Exception:  # noqa: BLE001
                pass
        if os.path.exists(local) and os.path.getsize(local) > 0:
            frames.append(_parse_product_zip(open(local, "rb").read()))
    # Recent (last ~1.5 yr) — predictable name, may not exist for closed stations.
    akt = f"tageswerte_KL_{station_id}_akt.zip"
    local_akt = os.path.join(cache_dir, akt)
    if not (os.path.exists(local_akt) and os.path.getsize(local_akt) > 0):
        try:
            rr = requests.get(_DWD_KL + "recent/" + akt, timeout=120)
            if rr.ok and rr.content[:2] == b"PK":
                open(local_akt, "wb").write(rr.content)
        except Exception:  # noqa: BLE001
            pass
    if os.path.exists(local_akt) and os.path.getsize(local_akt) > 0:
        frames.append(_parse_product_zip(open(local_akt, "rb").read()))

    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True).drop_duplicates("DATE").sort_values("DATE")
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# .WTH writer (TESTABLE with synthetic data; no network)
# ---------------------------------------------------------------------------

def _write_wth(df: pd.DataFrame, pid: str, lat: float, lon: float,
               elev, output_dir: str) -> str:
    """Write one DSSAT .WTH from a daily DataFrame.

    *df* must contain DATE (YYYYDOY str), SRAD, TMAX, TMIN, RAIN, TDEW, RH2M, WIND.
    """
    tav = _calc_tav(df)
    amp = _calc_amp(df)
    elev_str = f"{elev:5.0f}" if (elev is not None and np.isfinite(elev)) else "  -99"
    header = (
        f"$WEATHER DATA: DWD (Point ID: {pid})\n"
        f"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
        f"  DWD  {lat:8.4f} {lon:8.4f} {elev_str} {tav:5.1f} {amp:5.1f}   2.0  10.0\n"
        f"@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND"
    )
    lines = []
    for _, row in df.iterrows():
        line = (
            f"{row['DATE']:>7s}"
            f"{row['SRAD']:6.1f}{row['TMAX']:6.1f}{row['TMIN']:6.1f}"
            f"{row['RAIN']:6.1f}{row['TDEW']:6.1f}{row['RH2M']:6.1f}"
            f"{row['WIND']:6.1f}"
        )
        line = line.replace(" -99.0", "   -99")
        lines.append(line)
    out_path = os.path.join(output_dir, f"{pid}.WTH")
    with open(out_path, "w") as fh:
        fh.write(header + "\n")
        fh.write("\n".join(lines) + "\n")
    return out_path


def _build_point_frame(station_daily: pd.DataFrame, lat: float,
                       start_year: int, end_year: int) -> pd.DataFrame:
    """Map a station's raw daily product to a DSSAT-ready frame for [start,end]."""
    d = station_daily.copy()
    d = d[(d["DATE"].dt.year >= start_year) & (d["DATE"].dt.year <= end_year)]
    if d.empty:
        return d
    doy = d["DATE"].dt.dayofyear.values
    out = pd.DataFrame({
        "DATE": [f"{t.year}{t.dayofyear:03d}" for t in d["DATE"]],
        "YEAR": d["DATE"].dt.year.values,
        "MM": d["DATE"].dt.month.values,
        "TMAX": pd.to_numeric(d.get("TXK"), errors="coerce").values,
        "TMIN": pd.to_numeric(d.get("TNK"), errors="coerce").values,
        "RAIN": pd.to_numeric(d.get("RSK"), errors="coerce").values,
        "RH2M": pd.to_numeric(d.get("UPM"), errors="coerce").values,
        "WIND": pd.to_numeric(d.get("FM"), errors="coerce").values,
    })
    sdk = pd.to_numeric(d.get("SDK"), errors="coerce").values
    out["SRAD"] = _srad_from_sunshine(lat, doy, sdk)
    out["TDEW"] = _tdew_from_vapour_pressure(pd.to_numeric(d.get("VPM"), errors="coerce").values)
    # Drop rows with no temperature at all (station gaps); keep partial-var days.
    out = out[out["TMAX"].notna() & out["TMIN"].notna()].reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# Public entry point (mirrors the other weather sources + dwd_cache_dir)
# ---------------------------------------------------------------------------

def process_weather_dwd(
    shapefile,           # GeoDataFrame
    start_year: int,
    end_year: int,
    output_dir: str,
    id_col: str,
    lat_col: str,
    lon_col: str,
    n_cores: int,        # kept for API compatibility (fetch is shared/cached)
    log_file: str,
    dwd_cache_dir: str,
    max_station_km: float = 70.0,
) -> None:
    """Build DSSAT .WTH files from the nearest DWD station for each grid point.

    For each point the nearest station whose record overlaps [start_year,end_year]
    (within *max_station_km*) is used. Station downloads are cached and shared, so
    many nearby points reuse one fetch. SRAD is estimated from sunshine duration.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(dwd_cache_dir, exist_ok=True)
    end_year = min(end_year, date.today().year)

    print(f"--- Starting DWD Download (Years: {start_year}–{end_year}) ---")
    stations = _dwd_stations(dwd_cache_dir)
    hist_idx = _historical_index(dwd_cache_dir)
    # Keep only stations whose record overlaps the requested window.
    stations = stations[(stations["von"] // 10000 <= end_year)
                        & (stations["bis"] // 10000 >= start_year)].reset_index(drop=True)
    if stations.empty:
        print("DWD: no stations cover the requested period."); return
    print(f"DWD: {len(stations)} candidate stations cover {start_year}-{end_year}.")

    pts = shapefile.copy()
    if hasattr(pts, "geometry"):
        pts = pts.to_crs("EPSG:4326")
        pts[lat_col] = pts.geometry.y
        pts[lon_col] = pts.geometry.x

    st_lat = stations["lat"].values
    st_lon = stations["lon"].values

    def _nearest_station(lat, lon):
        # Great-circle distance to every candidate station (km).
        dlat = np.radians(st_lat - lat)
        dlon = np.radians(st_lon - lon)
        a = (np.sin(dlat / 2) ** 2 + np.cos(np.radians(lat))
             * np.cos(np.radians(st_lat)) * np.sin(dlon / 2) ** 2)
        dist = 6371.0 * 2 * np.arcsin(np.sqrt(a))
        order = np.argsort(dist)
        return order, dist

    station_cache = {}   # station_id -> parsed daily DataFrame (shared across points)

    def _get_station_daily(sid):
        if sid not in station_cache:
            station_cache[sid] = _fetch_station(sid, hist_idx.get(sid), dwd_cache_dir)
        return station_cache[sid]

    written = 0
    for _, prow in pts.iterrows():
        pid = str(prow[id_col])
        lat = float(prow[lat_col]); lon = float(prow[lon_col])
        out_path = os.path.join(output_dir, f"{pid}.WTH")
        if os.path.exists(out_path):
            written += 1
            continue
        try:
            order, dist = _nearest_station(lat, lon)
            frame = None
            used = None
            for k in order[:8]:                       # try up to 8 nearest stations
                if dist[k] > max_station_km:
                    break
                sid = stations.iloc[k]["station_id"]
                daily = _get_station_daily(sid)
                if daily.empty:
                    continue
                f = _build_point_frame(daily, stations.iloc[k]["lat"], start_year, end_year)
                if not f.empty and len(f) >= 30:      # need a usable record
                    frame = f; used = stations.iloc[k]; break
            if frame is None:
                raise ValueError(f"no DWD station within {max_station_km:.0f} km with data for {start_year}-{end_year}")
            frame = frame.fillna(-99)
            elev = used["elev"] if "elev" in used else np.nan
            # Write at the STATION coordinates (where the obs actually are).
            _write_wth(frame, pid, float(used["lat"]), float(used["lon"]), elev, output_dir)
            written += 1
        except Exception as exc:  # noqa: BLE001
            msg = f"\n--- ERROR ---\nDWD point {pid} ({lat:.3f},{lon:.3f}): {exc}\n"
            print(msg)
            with open(log_file, "a") as lf:
                lf.write(msg)

    print(f"\nDWD processing complete: {written}/{len(pts)} points written to '{output_dir}'.\n")
