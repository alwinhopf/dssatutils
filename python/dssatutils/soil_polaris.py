# File: soil_polaris.py
# Python port of soil_polaris.R
#
# Fetches POLARIS soil-property data (30 m probabilistic disaggregation of
# SSURGO over CONUS) via GeoTIFF tiles streamed with GDAL /vsicurl, derives
# DSSAT soil physics from POLARIS's published van Genuchten retention curve,
# and writes individual per-point .SOL files plus a mapping CSV.
#
# This is the Tier-0 deterministic drop-in: one profile per point built from a
# single statistic layer (default the median, ``p50``). The percentile layers
# (p5/p95) for soil-input uncertainty are intentionally out of scope here -- the
# ``stat`` argument is the seam a later ensemble/Monte-Carlo layer hooks into.
#
# References:
#   Chaney et al. (2016) Geoderma 274:54-67   (POLARIS 30 m soil series)
#   Chaney et al. (2019) Water Resour. Res. 55:2916-2938  (POLARIS soil props)
# Data root (keyless): http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/

import math
import os
import warnings
from datetime import datetime
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import requests

# --------------------------------------------------------------------------- #
# 0. Product layout & units
# --------------------------------------------------------------------------- #
POLARIS_BASE = "http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0"

# Properties needed to build a DSSAT profile. POLARIS publishes more (lambda,
# hb -- Brooks-Corey); we use the van Genuchten set + texture/bd/om/ph/ksat.
_VARS = ["clay", "sand", "silt", "bd", "om", "ph",
         "theta_r", "theta_s", "alpha", "n", "ksat"]

# These layers are stored as log10(value) and must be back-transformed (10**x).
#   om    -> log10(%)            alpha -> log10(kPa^-1)
#   ksat  -> log10(cm/hr)        hb    -> log10(kPa)   (unused here)
_LOG10_VARS = {"om", "alpha", "ksat", "hb"}

# POLARIS depth intervals (GlobalSoilMap standard), bottoms in cm.
_DEPTH_LABELS = ["0_5", "5_15", "15_30", "30_60", "60_100", "100_200"]
_DEPTH_BOTTOMS = [5, 15, 30, 60, 100, 200]
_DEPTH_CENTERS = [2.5, 10.0, 22.5, 45.0, 80.0, 150.0]

# Matric potentials (kPa, magnitudes) at which DUL / LL are evaluated on the
# van Genuchten curve. alpha is in kPa^-1 so psi must be in kPa to match.
# NOTE: if a future POLARIS release changes alpha's units to cm^-1, convert
# these to cm head (1 kPa ~ 10.2 cm) -- this is the #1 units gotcha.
PSI_DUL_KPA = 33.0      # field capacity (1/3 bar)
PSI_LL_KPA = 1500.0     # permanent wilting point (15 bar)

# Crash guards (same physical floors used by soil_soilgrids_online / soil_ssurgo
# to stop DSSAT's water balance SIGFPE-ing on degenerate sandy layers).
_LL_FLOOR = 0.02
_PAW_MIN = 0.04         # minimum DUL-LL and SSAT-DUL gap


# --------------------------------------------------------------------------- #
# 1. Tile addressing & value transforms (pure -> unit-testable offline)
# --------------------------------------------------------------------------- #
def _polaris_tile(lat: float, lon: float) -> str:
    """Return the 1 deg x 1 deg POLARIS tile token covering (lat, lon).

    Tiles are named ``lat{S}{N}_lon{W}{E}`` with integer bounds, e.g. a point at
    (42.35, -93.40) -> ``lat4243_lon-94-93``; (35.91, -101.40) ->
    ``lat3536_lon-102-101``.
    """
    s = math.floor(lat)
    w = math.floor(lon)
    return f"lat{s}{s + 1}_lon{w}{w + 1}"


def _backtransform(var: str, value: float) -> float:
    """Undo POLARIS log10 storage for the log-scaled properties."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return float("nan")
    return float(10.0 ** value) if var in _LOG10_VARS else float(value)


def _vg_theta(psi_kpa: float, theta_r: float, theta_s: float,
              alpha: float, n: float) -> float:
    """van Genuchten (1980) water content at matric potential *psi_kpa*.

    theta(psi) = theta_r + (theta_s - theta_r) / [1 + (alpha*|psi|)^n]^(1-1/n)
    *alpha* is already back-transformed to kPa^-1. Returns NaN if params are
    unusable (n<=1 makes m<=0).
    """
    if not (n and n > 1.0) or any(pd.isna(v) for v in (theta_r, theta_s, alpha)):
        return float("nan")
    m = 1.0 - 1.0 / n
    return theta_r + (theta_s - theta_r) / (1.0 + (alpha * abs(psi_kpa)) ** n) ** m


def _saxton_rawls(sand_pct: float, clay_pct: float, om_pct: float) -> Dict[str, float]:
    """Saxton & Rawls (2006) fallback when POLARIS van Genuchten params are
    missing for a layer. Identical formulation to soil_soilgrids_online.py."""
    S, C, OM = sand_pct / 100.0, clay_pct / 100.0, om_pct / 100.0
    t1500 = (-0.024 * S + 0.487 * C + 0.006 * OM + 0.005 * S * OM
             - 0.013 * C * OM + 0.068 * S * C + 0.031)
    lll = t1500 + (0.14 * t1500 - 0.02)
    t33 = (-0.251 * S + 0.195 * C + 0.011 * OM + 0.006 * S * OM
           - 0.027 * C * OM + 0.452 * S * C + 0.299)
    dul = t33 + (1.283 * t33 ** 2 - 0.374 * t33 - 0.015)
    ts33t = (0.278 * S + 0.034 * C + 0.022 * OM - 0.018 * S * OM
             - 0.027 * C * OM - 0.584 * S * C + 0.078)
    ts33 = ts33t + (0.636 * ts33t - 0.107)
    sat = dul + ts33 - 0.097 * S + 0.043
    return {"SLLL": lll, "SDUL": dul, "SSAT": sat}


def water_limits(theta_r, theta_s, alpha, n, *,
                 sand=None, clay=None, om_pct=None) -> Dict[str, float]:
    """DSSAT SLLL / SDUL / SSAT from the POLARIS van Genuchten curve, with a
    Saxton-Rawls fallback and the standard ordering/PAW crash guards applied.

    SLLL = theta(1500 kPa), SDUL = theta(33 kPa), SSAT = theta_s. The guards
    enforce SLLL >= 0.02 and SLLL < SDUL < SSAT with a >= 0.04 gap, which DSSAT
    requires (else IPSOIL / water-balance crashes -- see README C/B fixes).
    """
    lll = _vg_theta(PSI_LL_KPA, theta_r, theta_s, alpha, n)
    dul = _vg_theta(PSI_DUL_KPA, theta_r, theta_s, alpha, n)
    sat = float(theta_s) if not pd.isna(theta_s) else float("nan")

    if any(pd.isna(v) for v in (lll, dul, sat)):
        if None not in (sand, clay, om_pct) and not any(
                pd.isna(v) for v in (sand, clay, om_pct)):
            sr = _saxton_rawls(float(sand), float(clay), float(om_pct))
            lll = sr["SLLL"] if pd.isna(lll) else lll
            dul = sr["SDUL"] if pd.isna(dul) else dul
            sat = sr["SSAT"] if pd.isna(sat) else sat
        else:
            raise ValueError("No usable van Genuchten or texture data for layer.")

    lll = max(float(lll), _LL_FLOOR)
    dul = max(float(dul), lll + _PAW_MIN)
    sat = max(float(sat), dul + _PAW_MIN)
    return {"SLLL": lll, "SDUL": dul, "SSAT": sat}


def _ssks_cmhr(ksat: float) -> float:
    """DSSAT SSKS (cm/hr), clamped to the column range. *ksat* is already
    back-transformed to linear cm/hr by ``_backtransform``."""
    if pd.isna(ksat):
        return -99.0
    return float(min(999.0, max(0.0, ksat)))


# --------------------------------------------------------------------------- #
# 2. DSSAT .SOL formatter (mirrors soil_soilgrids_online._format_dssat_sol_file)
# --------------------------------------------------------------------------- #
def _format_dssat_sol_file(site_data: pd.DataFrame, output_dir: str,
                           source_tag: str = "p50",
                           source_name: str = "POLARIS v1.0") -> None:
    """Write one DSSAT .SOL file from per-layer POLARIS-derived properties."""
    if site_data.empty:
        raise ValueError("No soil layers found for this ID.")
    # Columns actually written to the profile; any all-NA among these is fatal.
    critical = ["clay", "silt", "bd"]
    if site_data[critical].isna().all().any():
        raise ValueError("Critical soil data (clay/silt/bulk density) all NA.")

    soil_id = str(site_data["ID"].iloc[0])
    lat = float(site_data["latitude"].iloc[0])
    lon = float(site_data["longitude"].iloc[0])
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"{soil_id}.SOL")

    lines = [
        f"*SOILS: {source_name}",
        f"! Source: {source_name} (statistic={source_tag})",
        "",
        f"*{soil_id[:10]:<10s}  POLARIS       {lat:9.3f} {lon:9.3f}",
        "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
        f" {soil_id[:11]:<11s} USA           {lat:9.3f} {lon:9.3f} ",
        "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
        "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001",
        "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC",
    ]

    site_data = site_data.sort_values("depth_bottom").reset_index(drop=True)
    for _, layer in site_data.iterrows():
        dc = float(layer["depth_center"])
        srgf = max(0.0, math.exp(-0.02 * dc))
        if srgf < 0.02:
            srgf = 0.0
        slhw = layer["ph"] if not pd.isna(layer["ph"]) else -99.0
        lines.append(
            f"{int(layer['depth_bottom']):6d}   -99"
            f" {layer['SLLL']:5.3f} {layer['SDUL']:5.3f} {layer['SSAT']:5.3f}"
            f" {srgf:5.2f} {layer['SSKS']:5.1f} {layer['bd']:5.2f} {layer['oc_pct']:5.2f}"
            f" {layer['clay']:5.1f} {layer['silt']:5.1f}   -99"
            f"   -99 {slhw:5.1f}   -99   -99   -99"
        )
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
# 3. Tile fetch (GDAL /vsicurl streaming, optional local cache)
# --------------------------------------------------------------------------- #
def _tile_source(var: str, stat: str, depth_label: str, tile: str,
                 cache_dir: Optional[str]) -> str:
    """Return a rasterio-openable path for one (var, stat, depth, tile) GeoTIFF.

    With *cache_dir* the tile is downloaded once and reused (HPC-friendly);
    otherwise it is streamed in place via GDAL /vsicurl (block-level range reads,
    so a small field never downloads the whole ~50 MB tile).
    """
    url = f"{POLARIS_BASE}/{var}/{stat}/{depth_label}/{tile}.tif"
    if not cache_dir:
        return f"/vsicurl/{url}"
    local = os.path.join(cache_dir, var, stat, depth_label, f"{tile}.tif")
    if not os.path.exists(local):
        os.makedirs(os.path.dirname(local), exist_ok=True)
        with requests.get(url, stream=True, timeout=180) as r:
            r.raise_for_status()
            tmp = local + ".tmp"
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
            os.replace(tmp, local)
    return local


def _fetch_polaris(gridfile, id_col: str, stat: str,
                   cache_dir: Optional[str]) -> pd.DataFrame:
    """Sample all POLARIS properties at every grid point, grouped by tile.

    Returns a long DataFrame: ID, prop, depth_bottom, depth_center, value
    (log-scaled props already back-transformed to linear units).
    """
    try:
        import rasterio  # noqa: F401
    except ImportError as exc:
        raise ImportError("rasterio is required for POLARIS. "
                          "Install with: pip install rasterio") from exc
    import rasterio

    gdf = gridfile.to_crs("EPSG:4326")
    lons = gdf.geometry.x.values
    lats = gdf.geometry.y.values
    ids = list(gridfile[id_col])

    # group point indices by the tile that contains them
    tiles: Dict[str, List[int]] = {}
    for i, (la, lo) in enumerate(zip(lats, lons)):
        tiles.setdefault(_polaris_tile(la, lo), []).append(i)

    records = []
    for tile, idxs in tiles.items():
        pts = [(float(lons[i]), float(lats[i])) for i in idxs]
        for var in _VARS:
            for d, dlabel in enumerate(_DEPTH_LABELS):
                src_path = _tile_source(var, stat, dlabel, tile, cache_dir)
                try:
                    with rasterio.open(src_path) as src:
                        sampled = list(src.sample(pts, masked=True))
                except Exception as exc:  # missing tile (e.g. offshore) etc.
                    warnings.warn(f"POLARIS skip {var}/{stat}/{dlabel}/{tile}: {exc}")
                    continue
                for k, cell in zip(idxs, sampled):
                    v = cell[0] if cell is not None else None
                    if v is None or np.ma.is_masked(v):
                        val = float("nan")
                    else:
                        val = _backtransform(var, float(v))
                    records.append({
                        "ID": ids[k], "prop": var,
                        "depth_bottom": _DEPTH_BOTTOMS[d],
                        "depth_center": _DEPTH_CENTERS[d],
                        "value": val,
                    })
    return pd.DataFrame(records)


# --------------------------------------------------------------------------- #
# 4. Public entry point
# --------------------------------------------------------------------------- #
def process_soils_polaris(
    gridfile,                 # GeoDataFrame with .geometry and an id column
    soilfile_csv_path: str,
    output_sol_dir: str,
    id_col: str,
    stat: str = "p50",
    cache_dir: Optional[str] = None,
) -> None:
    """Fetch POLARIS soil properties, derive DSSAT physics, and write per-point
    .SOL files plus a mapping CSV. Mirrors ``process_soils_polaris`` in R.

    Parameters
    ----------
    gridfile : GeoDataFrame
        Grid points with geometry and an ``id_col`` column.
    soilfile_csv_path : str
        Output path for the ID -> SOIL_ID mapping CSV.
    output_sol_dir : str
        Directory to receive ``<ID>.SOL`` files.
    id_col : str
        Grid-point ID column.
    stat : str
        POLARIS statistic layer to use. Default ``"p50"`` (median) -- the
        Tier-0 deterministic drop-in. ``"mean"``, ``"p5"`` and ``"p95"`` are
        valid too; building an uncertainty ensemble from p5/p50/p95 is left to
        a higher layer.
    cache_dir : str, optional
        If given, POLARIS tiles are downloaded here once and reused.
    """
    if stat not in ("p50", "mean", "mode", "p5", "p95"):
        raise ValueError(f"Unknown POLARIS statistic: {stat!r}")
    if not hasattr(gridfile, "geometry"):
        raise TypeError("gridfile must be a GeoDataFrame (with .geometry column).")

    print(f"--- POLARIS extraction (statistic={stat}, CONUS 30 m) ---")
    grid_wgs84 = gridfile.to_crs("EPSG:4326").copy()
    grid_wgs84["lon_wgs84"] = grid_wgs84.geometry.x
    grid_wgs84["lat_wgs84"] = grid_wgs84.geometry.y

    long_df = _fetch_polaris(gridfile, id_col, stat, cache_dir)
    if long_df.empty:
        raise RuntimeError("No POLARIS data extracted (coords outside CONUS, or "
                           "the Duke server is unreachable).")

    wide = long_df.pivot_table(index=["ID", "depth_bottom", "depth_center"],
                               columns="prop", values="value",
                               aggfunc="first").reset_index()
    for v in _VARS:
        if v not in wide.columns:
            wide[v] = np.nan

    # drop points with no texture anywhere -> no usable physics
    input_ids = set(long_df["ID"].unique())
    usable = wide.groupby("ID")[["sand", "clay"]].transform(
        lambda s: s.notna().any()).any(axis=1)
    wide = wide[usable].reset_index(drop=True)
    skipped = sorted(input_ids - set(wide["ID"].unique()))
    if skipped:
        print(f"  [warn] {len(skipped)} point(s) had no POLARIS texture and were "
              f"skipped: {', '.join(map(str, skipped[:10]))}"
              + (" ..." if len(skipped) > 10 else ""))
    if wide.empty:
        raise RuntimeError("No usable POLARIS data for any point.")

    # derived organic carbon (om already back-transformed to % by _backtransform)
    wide["oc_pct"] = wide["om"] / 1.724
    wide["SSKS"] = wide["ksat"].apply(_ssks_cmhr)

    lim = wide.apply(lambda r: pd.Series(water_limits(
        r["theta_r"], r["theta_s"], r["alpha"], r["n"],
        sand=r["sand"], clay=r["clay"], om_pct=r["om"])), axis=1)
    wide = pd.concat([wide, lim], axis=1)

    coords = grid_wgs84[[id_col, "lon_wgs84", "lat_wgs84"]].rename(
        columns={id_col: "ID", "lon_wgs84": "longitude", "lat_wgs84": "latitude"})
    final_df = wide.merge(coords, on="ID", how="inner")

    # 1. mapping CSV (point ID == soil ID, like the other gridded sources)
    pd.DataFrame({id_col: grid_wgs84[id_col].values,
                  "SOIL_ID": grid_wgs84[id_col].values}).to_csv(
        soilfile_csv_path, index=False)

    # 2. per-point .SOL files with error logging
    os.makedirs(output_sol_dir, exist_ok=True)
    log_path = os.path.join(output_sol_dir, "soil_processing_errors.log")
    with open(log_path, "w") as lf:
        lf.write(f"Log started: {datetime.now()}\n")

    success = errors = 0
    for uid in final_df["ID"].unique():
        subset = final_df[final_df["ID"] == uid].copy()
        try:
            _format_dssat_sol_file(subset, output_sol_dir, source_tag=stat)
            success += 1
        except Exception as exc:
            print(f"SKIPPED: ID {uid} | {exc}")
            with open(log_path, "a") as lf:
                lf.write(f"ID: {uid} | Error: {exc}\n")
            errors += 1

    print(f"POLARIS processing complete. Success: {success}, Errors: {errors}")
    print(f"See {log_path} for any skipped points.")
