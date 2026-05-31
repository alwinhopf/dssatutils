# File: soil_ssurgo.py
# Python port of soil_ssurgo.R
#
# Queries USDA SSURGO via the Soil Data Access (SDA) REST API, calculates
# DSSAT-ready soil physics using Saxton & Rawls (2006), and writes individual
# per-point .SOL files plus a mapping CSV.
#
# SDA REST endpoint: https://sdmdataaccess.nrcs.usda.gov/Tabular/post.rest
# Spatial SQL function: SDA_Get_Mukey_from_intersection_with_WktWgs84

import math
import os
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import numpy as np
import pandas as pd
import requests

try:
    from tqdm import tqdm
    _HAS_TQDM = True
except ImportError:
    _HAS_TQDM = False

_SDA_URL = "https://sdmdataaccess.nrcs.usda.gov/Tabular/post.rest"

# Standard SSURGO layer depth ranges (cm top, cm bottom)
_LAYER_RANGES = [
    (0, 5), (5, 20), (20, 35), (35, 50), (50, 65), (65, 80),
    (80, 95), (95, 110), (110, 125), (125, 140), (140, 155),
    (155, 170), (170, 185), (185, 200),
]


# ---------------------------------------------------------------------------
# SDA query helpers
# ---------------------------------------------------------------------------

def _sda_query(sql: str, max_retries: int = 3, delay: float = 5.0) -> Optional[pd.DataFrame]:
    """POST a SQL query to SDA and return a DataFrame, or None on failure."""
    for attempt in range(max_retries):
        try:
            r = requests.post(
                _SDA_URL,
                data={"query": sql, "format": "JSON+OBJECTS"},
                timeout=120,
            )
            r.raise_for_status()
            payload = r.json()
            table = payload.get("Table")
            if not table:
                return None
            rows = table[1:]  # first row is headers
            cols = table[0]
            return pd.DataFrame(rows, columns=cols)
        except Exception as exc:
            if attempt < max_retries - 1:
                time.sleep(delay)
            else:
                warnings.warn(f"SDA query failed: {exc}")
                return None


def _sda_spatial_mukeys(lat: float, lon: float,
                         max_retries: int = 3, delay: float = 5.0) -> Optional[list]:
    """Return list of mukeys intersecting the given WGS84 point."""
    wkt = f"POINT({lon} {lat})"
    sql = f"SELECT mukey FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('{wkt}')"
    df = _sda_query(sql, max_retries=max_retries, delay=delay)
    if df is None or df.empty:
        return None
    return df["mukey"].dropna().tolist()


def _format_in(values) -> str:
    """Format a list as SQL IN clause: ('a','b','c')."""
    cleaned = list({str(v) for v in values if v is not None})
    if not cleaned:
        return "('')"
    inner = ",".join(f"'{v}'" for v in cleaned)
    return f"({inner})"


# ---------------------------------------------------------------------------
# Soil physics: Saxton & Rawls (2006)
# ---------------------------------------------------------------------------

def _saxton_rawls(sand_pct: float, clay_pct: float, om_pct: float
                  ) -> tuple[float, float, float]:
    """
    Compute SLLL (wilting point), SDUL (field capacity), SSAT (saturation)
    as volumetric water fractions using Saxton & Rawls (2006) equations.
    """
    S = sand_pct / 100.0
    C = clay_pct / 100.0
    OM = om_pct / 100.0

    # Theta at 1500 kPa (wilting point)
    theta_1500t = (-0.024 * S + 0.487 * C + 0.006 * OM
                   + 0.005 * S * OM - 0.013 * C * OM
                   + 0.068 * S * C + 0.031)
    SLLL = theta_1500t + (0.14 * theta_1500t - 0.02)

    # Theta at 33 kPa (field capacity)
    theta_33t = (-0.251 * S + 0.195 * C + 0.011 * OM
                 + 0.006 * S * OM - 0.027 * C * OM
                 + 0.452 * S * C + 0.299)
    SDUL = theta_33t + (1.283 * theta_33t**2 - 0.374 * theta_33t - 0.015)

    # Saturation
    theta_s33t = (0.278 * S + 0.034 * C + 0.022 * OM
                  - 0.018 * S * OM - 0.027 * C * OM
                  - 0.584 * S * C + 0.078)
    theta_s33 = theta_s33t + (0.636 * theta_s33t - 0.107)
    SSAT = SDUL + theta_s33 - 0.097 * S + 0.043

    return float(SLLL), float(SDUL), float(SSAT)


# ---------------------------------------------------------------------------
# Weighted soil property aggregation per layer
# ---------------------------------------------------------------------------

def _calc_layer_props(props_df: pd.DataFrame,
                      top_cm: float, bot_cm: float) -> Optional[pd.Series]:
    """
    Compute component-percentage-weighted average clay, sand, OM, and bulk
    density for the depth interval [top_cm, bot_cm].  Returns None if no
    horizons overlap.
    """
    df = props_df.copy()
    df["adj_top"] = df["hzdept_r"].clip(lower=top_cm)
    df["adj_bot"] = df["hzdepb_r"].clip(upper=bot_cm)
    df["thickness"] = (df["adj_bot"] - df["adj_top"]).clip(lower=0)
    df = df[df["thickness"] > 0].copy()
    if df.empty:
        return None

    w = df["thickness"] * df["comppct_r"]
    total_w = w.sum()
    if total_w == 0:
        return None

    return pd.Series({
        "clay_pct": (df["claytotal_r"] * w).sum() / total_w,
        "sand_pct": (df["sandtotal_r"] * w).sum() / total_w,
        "om_pct":   (df["om_r"] * w).sum() / total_w,
        "bulk_density": (df["dbthirdbar_r"] * w).sum() / total_w,
    })


# ---------------------------------------------------------------------------
# DSSAT .SOL file writer
# ---------------------------------------------------------------------------

def _write_sol(profile: pd.DataFrame, output_dir: str) -> None:
    """Write one DSSAT .SOL file from a profile DataFrame."""
    soil_id = str(profile["ID"].iloc[0])
    lat = profile["latitude"].iloc[0]
    lon = profile["longitude"].iloc[0]
    path = os.path.join(output_dir, f"{soil_id}.SOL")
    if os.path.exists(path):
        return

    lines = [
        "*SOILS: USA SSURGO Soil Profiles",
        "! Generated from SSURGO database",
        "",
        f"*{soil_id:<6s}  SSURGO        {lat:9.3f} {lon:9.3f}",
        "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
        f" {soil_id:<11s} USA         {lat:9.3f} {lon:9.3f} ",
        "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
        "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001",
        "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC",
    ]

    for _, layer in profile.sort_values("depth_bottom").iterrows():
        slll = f"{layer['SLLL']:5.3f}".lstrip("0") or "0.000"
        sdul = f"{layer['SDUL']:5.3f}".lstrip("0") or "0.000"
        ssat = f"{layer['SSAT']:5.3f}".lstrip("0") or "0.000"
        # Ensure leading space format like R's sub("^0", " ", ...)
        slll = " " + slll if slll[0].isdigit() else slll
        sdul = " " + sdul if sdul[0].isdigit() else sdul
        ssat = " " + ssat if ssat[0].isdigit() else ssat

        depth = int(layer["depth_bottom"])
        depth_str = f"{depth:6d}" if depth >= 10 else f"{depth:6d}"
        om_sloc = layer["om_pct"] / 1.724  # OM → SOC
        lines.append(
            f"{depth_str}   -99 {slll} {sdul} {ssat}  1.00   -99"
            f" {layer['bulk_density']:5.2f} {om_sloc:5.2f}"
            f" {layer['clay_pct']:5.1f} {layer['silt_pct']:5.1f}"
            f"   -99   -99   -99   -99   -99   -99"
        )
    lines.append("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Per-point worker
# ---------------------------------------------------------------------------

def _process_point(args: dict) -> Optional[pd.DataFrame]:
    """Query SSURGO for one point and write its .SOL file."""
    ID = args["ID"]
    lat = args["lat"]
    lon = args["lon"]
    output_dir = args["output_dir"]

    # Skip if already written
    if os.path.exists(os.path.join(output_dir, f"{ID}.SOL")):
        return None

    # 1. Spatial query → mukeys
    mukeys = _sda_spatial_mukeys(lat, lon)
    if not mukeys:
        return None

    # 2. Bedrock depth
    q_bed = (
        f"SELECT mukey, brockdepmin FROM muaggatt "
        f"WHERE mukey IN {_format_in(mukeys)}"
    )
    bed_df = _sda_query(q_bed)
    bedrock_depth = 200.0
    if bed_df is not None and not bed_df.empty:
        brd = pd.to_numeric(bed_df["brockdepmin"], errors="coerce").dropna()
        if not brd.empty:
            bedrock_depth = float(brd.min())
    if not math.isfinite(bedrock_depth):
        bedrock_depth = 200.0

    # 3. Determine valid layers up to bedrock
    valid_layers = [(t, min(b, bedrock_depth))
                    for t, b in _LAYER_RANGES if t < bedrock_depth]
    if not valid_layers:
        valid_layers = [(0.0, bedrock_depth)]

    # 4. Query horizon properties
    q_prop = (
        "SELECT component.mukey, component.cokey, component.comppct_r, "
        "chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.claytotal_r, "
        "chorizon.sandtotal_r, chorizon.om_r, chorizon.dbthirdbar_r "
        "FROM component INNER JOIN chorizon ON component.cokey = chorizon.cokey "
        f"WHERE component.mukey IN {_format_in(mukeys)}"
    )
    props_df = _sda_query(q_prop)
    if props_df is None or props_df.empty:
        return None

    # Coerce numeric columns
    for col in ["hzdept_r", "hzdepb_r", "claytotal_r", "sandtotal_r",
                "om_r", "dbthirdbar_r", "comppct_r"]:
        props_df[col] = pd.to_numeric(props_df[col], errors="coerce")

    # 5. Aggregate per layer
    layer_rows = []
    for (top, bot) in valid_layers:
        agg = _calc_layer_props(props_df, top, bot)
        if agg is None:
            continue
        clay = float(agg["clay_pct"]) if not np.isnan(agg["clay_pct"]) else 20.0
        sand = float(agg["sand_pct"]) if not np.isnan(agg["sand_pct"]) else 40.0
        om   = float(agg["om_pct"])   if not np.isnan(agg["om_pct"])   else 1.0
        bd   = float(agg["bulk_density"]) if not np.isnan(agg["bulk_density"]) else 1.4
        silt = max(0.0, 100.0 - clay - sand)

        SLLL, SDUL, SSAT = _saxton_rawls(sand, clay, om)
        layer_rows.append({
            "ID": ID, "latitude": lat, "longitude": lon,
            "depth_top": top, "depth_bottom": bot,
            "clay_pct": clay, "sand_pct": sand, "silt_pct": silt,
            "om_pct": om, "bulk_density": bd,
            "SLLL": SLLL, "SDUL": SDUL, "SSAT": SSAT,
        })

    if not layer_rows:
        return None

    profile_df = pd.DataFrame(layer_rows)
    _write_sol(profile_df, output_dir)
    return profile_df


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def process_soils_ssurgo(
    grid_points,              # GeoDataFrame
    output_dir_csv: str,
    output_dir_individual: str,
    n_cores: int,
    id_col: str,
    lat_col: str,
    long_col: str,
    format_sql_func=None,     # kept for API compat; unused internally
) -> bool:
    """
    Query SSURGO for every point in *grid_points*, write per-point .SOL files
    to *output_dir_individual*, and append a mapping CSV to *output_dir_csv*.
    Mirrors the R ``process_soils_ssurgo`` function including smart-resume.
    """
    import geopandas as gpd

    print("Starting SSURGO Processing (Smart Resume Mode)...")
    os.makedirs(output_dir_individual, exist_ok=True)

    # --- Smart resume: skip points that already have a .SOL file ---
    existing = {
        os.path.splitext(f)[0]
        for f in os.listdir(output_dir_individual)
        if f.endswith(".SOL")
    }

    # Build task list from GeoDataFrame
    gdf = grid_points.copy()
    if hasattr(gdf, "geometry"):
        import geopandas as gpd
        gdf = gdf.to_crs("EPSG:4326")
        gdf[lat_col] = gdf.geometry.y
        gdf[long_col] = gdf.geometry.x

    all_ids = gdf[id_col].astype(str).tolist()
    missing_mask = [str(pid) not in existing for pid in all_ids]
    to_process = gdf[missing_mask].reset_index(drop=True)

    n_total = len(all_ids)
    n_skip = n_total - len(to_process)
    n_proc = len(to_process)
    print(f"Resume Check: Found {n_skip} existing profiles. Processing {n_proc} remaining.")

    if n_proc == 0:
        print("All soil profiles already exist. Skipping SSURGO processing.")
        return True

    CHUNK_SIZE = 10_000
    num_chunks = math.ceil(n_proc / CHUNK_SIZE)
    print(f"Processing {n_proc} points in {num_chunks} chunk(s)...")

    csv_header_written = os.path.exists(output_dir_csv)

    for chunk_i in range(num_chunks):
        s = chunk_i * CHUNK_SIZE
        e = min((chunk_i + 1) * CHUNK_SIZE, n_proc)
        chunk = to_process.iloc[s:e]
        print(f"  > Chunk {chunk_i+1}/{num_chunks} (Points {s+1} – {e})")

        tasks = [
            {"ID": str(row[id_col]),
             "lat": float(row[lat_col]),
             "lon": float(row[long_col]),
             "output_dir": output_dir_individual}
            for _, row in chunk.iterrows()
        ]

        results = []
        iter_obj = (
            tqdm(tasks, desc=f"Chunk {chunk_i+1}", unit="pt") if _HAS_TQDM else tasks
        )

        with ThreadPoolExecutor(max_workers=min(n_cores, 16)) as pool:
            future_map = {pool.submit(_process_point, t): t["ID"] for t in tasks}
            for fut in as_completed(future_map):
                pid = future_map[fut]
                try:
                    res = fut.result()
                    if res is not None:
                        results.append(res)
                except Exception as exc:
                    warnings.warn(f"Point {pid} failed: {exc}")

        if results:
            chunk_df = pd.concat(results, ignore_index=True)
            chunk_df.to_csv(
                output_dir_csv,
                mode="a",
                index=False,
                header=not csv_header_written,
            )
            csv_header_written = True

    print("SSURGO Processing Complete.")
    return True
