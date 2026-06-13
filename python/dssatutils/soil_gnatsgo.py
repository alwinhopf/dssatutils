# File: soil_gnatsgo.py
# Python port of soil_gnatsgo.R
#
# Queries USDA gNATSGO (gridded National Soil Survey Geographic Database) and
# writes per-point DSSAT .SOL files plus a mapping CSV.
#
# WHY gNATSGO over SSURGO: gNATSGO is a complete, gap-free 30 m soil grid for the
# conterminous US. It merges SSURGO + STATSGO2 + Raster Soil Surveys, so it has a
# map unit EVERYWHERE there is land — filling the "no-coverage" holes that plain
# SSURGO leaves (un-surveyed areas, some military/tribal land). For the vast
# majority of cropland the underlying map unit IS the SSURGO one, so values match
# SSURGO; the difference shows up only where SSURGO had a gap.
#
# ACCESS (two public services, no key required):
#   1. Map-unit key at a point — SoilWeb Web Coverage Service (UC Davis), the
#      same 30 m gNATSGO mukey grid the R `soilDB::mukey.wcs()` uses. Native CRS
#      EPSG:5070 (CONUS Albers).
#   2. Tabular soil properties for that mukey — USDA Soil Data Access (SDA), the
#      same REST endpoint soil_ssurgo.py uses (SDA hosts SSURGO + STATSGO2
#      tabular). A small fraction of gNATSGO mukeys are Raster-Soil-Survey keys
#      that are NOT in SDA; those are logged as "no-tabular".
#
# This module reuses the soil physics, layer aggregation, depth handling and
# failure-logging conventions of soil_ssurgo.py so SSURGO and gNATSGO profiles
# are byte-comparable where they share a map unit. Only the spatial lookup (WCS
# mukey grid instead of the SDA polygon intersect) differs.

import io
import math
import os
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests

from .soil_ssurgo import (
    _LAYER_RANGES,
    _calc_layer_props,
    _failure,
    _format_in,
    _mapunit_names,
    _saxton_rawls,
    _sda_query,
    _sda_query_result,
)

try:
    from tqdm import tqdm
    _HAS_TQDM = True
except ImportError:
    _HAS_TQDM = False

# SoilWeb WCS (UC Davis California Soil Resource Lab) — gNATSGO mukey grid.
# Same endpoint/mapfile/coverage the soilDB package targets.
_WCS_BASE = ("http://casoilresource.lawr.ucdavis.edu/cgi-bin/mapserv"
             "?map=/data1/website/wcs/mukey-grids.map")
_WCS_COVERAGE = "gnatsgo"
_WCS_CRS = "EPSG:5070"   # CONUS Albers, native gNATSGO grid CRS
_WCS_RES = 30.0          # native 30 m resolution


# ---------------------------------------------------------------------------
# Spatial lookup: gNATSGO mukey at a point via the SoilWeb WCS
# ---------------------------------------------------------------------------

def _gnatsgo_mukey(lat: float, lon: float, buffer_m: float = 45.0,
                   max_retries: int = 3, delay: float = 5.0) -> dict:
    """Return {ok, mukey, error} for the gNATSGO map-unit key at a WGS84 point.

    Fetches a tiny (a few-pixel) GeoTIFF window around the point from the WCS and
    samples the centre cell. mukey 0 / nodata means "no map unit here" (water,
    outside CONUS). Requires rasterio + pyproj.
    """
    try:
        import rasterio
        from pyproj import Transformer
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "mukey": None, "error": f"dependency missing (rasterio/pyproj): {exc}"}

    try:
        x, y = Transformer.from_crs(4326, 5070, always_xy=True).transform(lon, lat)
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "mukey": None, "error": f"reprojection failed: {exc}"}

    params = {
        "SERVICE": "WCS", "VERSION": "2.0.1", "REQUEST": "GetCoverage",
        "COVERAGEID": _WCS_COVERAGE, "FORMAT": "image/tiff",
        "SUBSETTINGCRS": _WCS_CRS,
        "SUBSET": [f"x({x - buffer_m},{x + buffer_m})",
                   f"y({y - buffer_m},{y + buffer_m})"],
        "RESOLUTION": [f"x({_WCS_RES})", f"y({_WCS_RES})"],
        "GEOTIFF:COMPRESSION": "DEFLATE",
    }
    last_error = None
    for attempt in range(max_retries):
        try:
            r = requests.get(_WCS_BASE, params=params, timeout=120)
            ct = r.headers.get("Content-Type", "")
            if "tiff" not in ct and r.content[:2] not in (b"II", b"MM"):
                # MapServer returns an HTML error page on a bad request.
                last_error = f"WCS returned non-raster ({ct}): {r.content[:160]!r}"
            else:
                with rasterio.open(io.BytesIO(r.content)) as ds:
                    arr = ds.read(1)
                if arr.size == 0:
                    return {"ok": True, "mukey": None, "error": None}
                cr, cc = arr.shape[0] // 2, arr.shape[1] // 2
                mukey = int(arr[cr, cc])
                nodata = {0, -2147483648, 2147483647}
                if mukey in nodata or mukey < 0:
                    return {"ok": True, "mukey": None, "error": None}
                return {"ok": True, "mukey": str(mukey), "error": None}
        except Exception as exc:  # noqa: BLE001
            last_error = str(exc)
        if attempt < max_retries - 1:
            import time
            time.sleep(delay)
    return {"ok": False, "mukey": None, "error": last_error or "unknown WCS error"}


# ---------------------------------------------------------------------------
# DSSAT .SOL file writer (gNATSGO-labelled; same column layout as SSURGO)
# ---------------------------------------------------------------------------

def _write_sol(profile: pd.DataFrame, output_dir: str) -> None:
    """Write one DSSAT .SOL file from a profile DataFrame (gNATSGO header)."""
    soil_id = str(profile["ID"].iloc[0])
    lat = profile["latitude"].iloc[0]
    lon = profile["longitude"].iloc[0]
    path = os.path.join(output_dir, f"{soil_id}.SOL")
    if os.path.exists(path):
        return

    lines = [
        "*SOILS: USA gNATSGO Soil Profiles",
        "! Generated from gNATSGO (SoilWeb WCS mukey grid + USDA SDA tabular)",
        "",
        f"*{soil_id:<6s}  gNATSGO       {lat:9.3f} {lon:9.3f}",
        "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
        f" {soil_id:<11s} USA         {lat:9.3f} {lon:9.3f} ",
        "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
        "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001",
        "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC",
    ]

    for _, layer in profile.sort_values("depth_bottom").iterrows():
        def _f3(v):
            s = f"{v:5.3f}"
            return (" " + s[1:]) if s.startswith("0.") else s
        slll, sdul, ssat = _f3(layer["SLLL"]), _f3(layer["SDUL"]), _f3(layer["SSAT"])
        depth = int(layer["depth_bottom"])
        om_sloc = layer["om_pct"] / 1.724  # OM → SOC
        lines.append(
            f"{depth:6d}   -99 {slll} {sdul} {ssat}  1.00   -99"
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

def _process_point(args: dict):
    """Look up the gNATSGO mukey for one point and build its .SOL file."""
    ID = args["ID"]
    lat = args["lat"]
    lon = args["lon"]
    output_dir = args["output_dir"]

    if os.path.exists(os.path.join(output_dir, f"{ID}.SOL")):
        return None

    # 1. Spatial lookup → gNATSGO mukey (WCS grid).
    res = _gnatsgo_mukey(lat, lon)
    if not res["ok"]:
        return _failure(ID, lat, lon, f"network: gNATSGO WCS lookup failed ({res['error']})")
    mukey = res["mukey"]
    if mukey is None:
        return _failure(ID, lat, lon, "no-coverage: no gNATSGO map unit at this location (water / outside CONUS)")
    mukeys = [mukey]

    # 2. Bedrock depth (SDA) — identical to SSURGO.
    bed_df = _sda_query(f"SELECT mukey, brockdepmin FROM muaggatt WHERE mukey IN {_format_in(mukeys)}")
    bedrock_depth = 200.0
    if bed_df is not None and not bed_df.empty:
        brd = pd.to_numeric(bed_df["brockdepmin"], errors="coerce").dropna()
        if not brd.empty:
            bedrock_depth = float(brd.min())
    if not math.isfinite(bedrock_depth):
        bedrock_depth = 200.0

    valid_layers = [(t, min(b, bedrock_depth)) for t, b in _LAYER_RANGES if t < bedrock_depth]
    if not valid_layers:
        valid_layers = [(0.0, bedrock_depth)]

    # 3. Horizon properties (SDA).
    q_prop = (
        "SELECT component.mukey, component.cokey, component.comppct_r, "
        "chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.claytotal_r, "
        "chorizon.sandtotal_r, chorizon.om_r, chorizon.dbthirdbar_r "
        "FROM component INNER JOIN chorizon ON component.cokey = chorizon.cokey "
        f"WHERE component.mukey IN {_format_in(mukeys)}"
    )
    prop_res = _sda_query_result(q_prop)
    if not prop_res["ok"]:
        return _failure(ID, lat, lon, f"network: SDA horizon query failed ({prop_res['error']})")
    props_df = prop_res["data"]
    if props_df is None or props_df.empty:
        # Either a non-soil map unit (Water/Urban/Rock) or a Raster-Soil-Survey
        # mukey that is not present in SDA. Distinguish via the mapunit table.
        names = _mapunit_names(mukeys)
        muname = "; ".join(sorted({v for v in names.values() if v and v.lower() != "nan"}))
        if not names:
            return _failure(ID, lat, lon,
                            f"no-tabular: gNATSGO mukey {mukey} not found in SDA (likely a Raster Soil Survey unit)")
        suffix = f" [{muname}]" if muname else ""
        return _failure(ID, lat, lon,
                        f"no-soil: map unit has no soil horizons{suffix} — typically Water / Urban / Pits / Rock outcrop (mukey {mukey})")

    for col in ["hzdept_r", "hzdepb_r", "claytotal_r", "sandtotal_r", "om_r", "dbthirdbar_r", "comppct_r"]:
        props_df[col] = pd.to_numeric(props_df[col], errors="coerce")

    # 4. Aggregate per layer + Saxton & Rawls (shared with SSURGO).
    layer_rows = []
    for (top, bot) in valid_layers:
        agg = _calc_layer_props(props_df, top, bot)
        if agg is None:
            continue
        clay = float(agg["clay_pct"]) if not np.isnan(agg["clay_pct"]) else 20.0
        sand = float(agg["sand_pct"]) if not np.isnan(agg["sand_pct"]) else 40.0
        om = float(agg["om_pct"]) if not np.isnan(agg["om_pct"]) else 1.0
        bd = float(agg["bulk_density"]) if not np.isnan(agg["bulk_density"]) else 1.4
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
        names = _mapunit_names(mukeys)
        muname = "; ".join(sorted({v for v in names.values() if v and v.lower() != "nan"}))
        return _failure(ID, lat, lon,
                        f"no-layers: horizon data present but no usable layers after depth filtering (bedrock {bedrock_depth:g} cm; muname {muname or 'unknown'})")

    profile_df = pd.DataFrame(layer_rows)
    _write_sol(profile_df, output_dir)
    return profile_df


# ---------------------------------------------------------------------------
# Public entry point — signature mirrors process_soils_ssurgo
# ---------------------------------------------------------------------------

def process_soils_gnatsgo(
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
    Resolve gNATSGO soil for every point in *grid_points*, write per-point .SOL
    files to *output_dir_individual*, and append a mapping CSV to *output_dir_csv*.
    Smart-resume: points that already have a .SOL are skipped. Mirrors
    ``process_soils_ssurgo`` exactly (same signature, chunking, failure log).
    """
    print("Starting gNATSGO Processing (Smart Resume Mode)...")
    os.makedirs(output_dir_individual, exist_ok=True)

    existing = {
        os.path.splitext(f)[0]
        for f in os.listdir(output_dir_individual)
        if f.endswith(".SOL")
    }

    gdf = grid_points.copy()
    if hasattr(gdf, "geometry"):
        gdf = gdf.to_crs("EPSG:4326")
        gdf[lat_col] = gdf.geometry.y
        gdf[long_col] = gdf.geometry.x

    all_ids = gdf[id_col].astype(str).tolist()
    missing_mask = [str(pid) not in existing for pid in all_ids]
    to_process = gdf[missing_mask].reset_index(drop=True)

    n_total = len(all_ids)
    n_proc = len(to_process)
    print(f"Resume Check: Found {n_total - n_proc} existing profiles. Processing {n_proc} remaining.")
    if n_proc == 0:
        print("All soil profiles already exist. Skipping gNATSGO processing.")
        return True

    CHUNK_SIZE = 10_000
    num_chunks = math.ceil(n_proc / CHUNK_SIZE)
    print(f"Processing {n_proc} points in {num_chunks} chunk(s)...")

    csv_header_written = os.path.exists(output_dir_csv)
    failures = []

    for chunk_i in range(num_chunks):
        s = chunk_i * CHUNK_SIZE
        e = min((chunk_i + 1) * CHUNK_SIZE, n_proc)
        chunk = to_process.iloc[s:e]
        print(f"  > Chunk {chunk_i + 1}/{num_chunks} (Points {s + 1} – {e})")

        tasks = [
            {"ID": str(row[id_col]),
             "lat": float(row[lat_col]),
             "lon": float(row[long_col]),
             "output_dir": output_dir_individual}
            for _, row in chunk.iterrows()
        ]

        results = []
        iter_obj = tqdm(tasks, desc=f"Chunk {chunk_i + 1}", unit="pt") if _HAS_TQDM else tasks

        # WCS likes modest concurrency; cap below the SSURGO ceiling.
        with ThreadPoolExecutor(max_workers=min(n_cores, 8)) as pool:
            future_map = {pool.submit(_process_point, t): t["ID"] for t in tasks}
            for fut in as_completed(future_map):
                pid = future_map[fut]
                try:
                    res = fut.result()
                    if isinstance(res, dict) and res.get("_fail"):
                        failures.append(res)
                    elif res is not None:
                        results.append(res)
                except Exception as exc:  # noqa: BLE001
                    failures.append(_failure(str(pid), math.nan, math.nan, f"network: worker failed ({exc})"))

        if results:
            chunk_df = pd.concat(results, ignore_index=True)
            chunk_df.to_csv(output_dir_csv, mode="a", index=False, header=not csv_header_written)
            csv_header_written = True

    if failures:
        fail_df = pd.DataFrame(failures)[["ID", "latitude", "longitude", "reason"]]
        base = os.path.splitext(os.path.basename(output_dir_csv))[0]
        failure_log = os.path.join(os.path.dirname(output_dir_csv), f"{base}_download_failures.csv")
        fail_df.to_csv(failure_log, index=False)

        categories = fail_df["reason"].str.split(":", n=1).str[0].fillna("unknown")
        counts = categories.value_counts().to_dict()
        labels = {
            "network": "WCS/SDA request failed after retries",
            "no-coverage": "no gNATSGO map unit here (water / outside CONUS)",
            "no-soil": "non-soil map unit (Water, Urban, Pits, Rock) — no horizons exist",
            "no-tabular": "gNATSGO mukey not in SDA (Raster Soil Survey unit)",
            "no-layers": "horizons present but unusable after depth filtering",
        }
        print(f"[gNATSGO] {len(fail_df)} of {n_proc} processed point(s) produced NO soil profile:")
        for key in sorted(counts):
            print(f"   - {key:<12} {counts[key]:4d}   ({labels.get(key, 'see failure CSV')})")
        print(f"   Per-point details (ID, lat, long, reason) -> {failure_log}")

    print("gNATSGO Processing Complete.")
    return True
