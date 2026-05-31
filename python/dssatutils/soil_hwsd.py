# File: soil_hwsd.py
# ---------------------------------------------------------------------------
# Soil source: FAO/IIASA Harmonized World Soil Database v2.0 (HWSD2).
#
# WHY: HWSD2 is the FAO "official" harmonized global soil database (~1 km) and
# the long-standing reference for GLOBAL gridded crop-model studies (e.g. the
# AgMIP/GGCMI intercomparisons historically used HWSD). It complements ISRIC
# SoilGrids with an independent, expert-harmonized product.
#
# ACCESS: HWSD2 is not a streaming API — download it once from FAO and point
# the pipeline at the two files (mirrors the SOILGRIDS_10K external-file model):
#   • HWSD2 raster of mapping-unit IDs (GeoTIFF / BIL): `hwsd_raster_file`
#   • HWSD2 attribute database (SQLite .sqlite/.db):     `hwsd_db_file`
#   FAO HWSD v2.0: https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v2-0/
#
# This samples the raster at each point to get the HWSD2 mapping-unit (SMU) ID,
# looks up the dominant soil component's layers in the SQLite DB, computes DSSAT
# physics (Saxton & Rawls 2006, shared with the SoilGrids module), and writes
# per-point .SOL files plus a mapping CSV. Points over no-data cells are skipped
# with a warning (never crash the run).
# ---------------------------------------------------------------------------

import os
import sqlite3
from datetime import datetime

import numpy as np
import pandas as pd

# Reuse the Saxton & Rawls physics and the DSSAT .SOL writer (single source).
from .soil_soilgrids_online import _calculate_soil_physics, _format_dssat_sol_file

# HWSD2 layer table is documented as HWSD2_LAYERS. Column names below are matched
# case-insensitively with fallbacks, so minor schema variants still work.
_LAYER_TABLE_CANDIDATES = ["HWSD2_LAYERS", "HWSD2_LAYER", "LAYERS", "D_LAYERS"]
_COL = {  # canonical -> list of accepted source column names (case-insensitive)
    "smu":    ["HWSD2_SMU_ID", "SMU_ID", "HWSD2_SMU", "SMU"],
    "seq":    ["SEQUENCE", "SEQ"],
    "share":  ["SHARE", "PERCENT", "PCT"],
    "top":    ["TOPDEP", "TOP_DEPTH", "TOP"],
    "bot":    ["BOTDEP", "BOT_DEPTH", "BOTTOM", "BOT"],
    "sand":   ["SAND", "SAND_PCT"],
    "silt":   ["SILT", "SILT_PCT"],
    "clay":   ["CLAY", "CLAY_PCT"],
    "bulk":   ["BULK", "BULK_DENSITY", "BD", "REF_BULK_DENSITY"],
    "oc":     ["ORG_CARBON", "OC", "ORGANIC_CARBON", "SOC"],
    "coarse": ["COARSE", "GRAVEL", "CFVO"],
}


def _resolve_columns(df: pd.DataFrame) -> dict:
    """Map canonical keys to the actual column names present in *df*."""
    lower = {c.lower(): c for c in df.columns}
    out = {}
    for key, candidates in _COL.items():
        for cand in candidates:
            if cand.lower() in lower:
                out[key] = lower[cand.lower()]
                break
    return out


def _find_layer_table(conn: sqlite3.Connection) -> str:
    names = pd.read_sql_query(
        "SELECT name FROM sqlite_master WHERE type IN ('table','view')", conn
    )["name"].tolist()
    lower = {n.lower(): n for n in names}
    for cand in _LAYER_TABLE_CANDIDATES:
        if cand.lower() in lower:
            return lower[cand.lower()]
    # Fall back to the first table that has a SAND-like column.
    for n in names:
        cols = pd.read_sql_query(f"SELECT * FROM '{n}' LIMIT 1", conn).columns
        if any(c.lower() == "sand" for c in cols):
            return n
    raise RuntimeError(
        f"Could not find an HWSD2 layer table in the DB. Tables: {names}")


def _sample_smu_ids(raster_file: str, lons, lats) -> np.ndarray:
    """Sample the HWSD2 mapping-unit raster at each (lon, lat). nodata -> -1."""
    try:
        import rasterio
        from pyproj import Transformer
    except ImportError as exc:
        raise ImportError(
            "rasterio and pyproj are required for HWSD. "
            "Install with: pip install rasterio pyproj") from exc

    with rasterio.open(raster_file) as src:
        # Reproject query points to the raster CRS if needed.
        dst_crs = src.crs.to_string() if src.crs else "EPSG:4326"
        if dst_crs not in ("EPSG:4326", "epsg:4326"):
            tr = Transformer.from_crs("EPSG:4326", dst_crs, always_xy=True)
            xs, ys = tr.transform(np.asarray(lons), np.asarray(lats))
        else:
            xs, ys = np.asarray(lons), np.asarray(lats)
        nodata = src.nodata
        out = []
        for v in src.sample(zip(xs, ys), masked=True):
            cell = v[0]
            if cell is np.ma.masked or np.ma.is_masked(cell):
                out.append(-1)
            elif nodata is not None and float(cell) == float(nodata):
                out.append(-1)
            else:
                out.append(int(cell))
        return np.array(out, dtype=int)


def process_soils_hwsd(
    grid_points,             # GeoDataFrame (or DataFrame with lon/lat columns)
    hwsd_raster_file: str,
    hwsd_db_file: str,
    output_csv_path: str,
    output_sol_dir: str,
    id_col: str = "ID",
    lat_col: str = "LAT",
    long_col: str = "LONG",
) -> None:
    """Build per-point DSSAT .SOL files from HWSD2 (raster + SQLite)."""
    print("--- Starting HWSD2 Extraction ---")
    for f, what in [(hwsd_raster_file, "raster"), (hwsd_db_file, "database")]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"HWSD2 {what} file not found: {f}")

    # Coordinates (support GeoDataFrame geometry or explicit lat/lon columns).
    if hasattr(grid_points, "geometry") and grid_points.geometry.notna().any():
        g = grid_points.to_crs("EPSG:4326") if hasattr(grid_points, "to_crs") else grid_points
        lons = g.geometry.x.values
        lats = g.geometry.y.values
    else:
        lons = grid_points[long_col].astype(float).values
        lats = grid_points[lat_col].astype(float).values
    ids = grid_points[id_col].astype(str).values

    # 1. Sample HWSD2 SMU IDs at each point.
    smu_ids = _sample_smu_ids(hwsd_raster_file, lons, lats)

    # 2. Load the layer table once; resolve schema columns.
    conn = sqlite3.connect(hwsd_db_file)
    try:
        table = _find_layer_table(conn)
        layers = pd.read_sql_query(f"SELECT * FROM '{table}'", conn)
    finally:
        conn.close()
    cols = _resolve_columns(layers)
    required = ["smu", "sand", "clay", "bot"]
    missing = [k for k in required if k not in cols]
    if missing:
        raise RuntimeError(
            f"HWSD2 layer table '{table}' is missing required columns "
            f"for: {missing}. Found columns: {list(layers.columns)}")

    # 3. Build per-point layer records (dominant component per SMU).
    os.makedirs(output_sol_dir, exist_ok=True)
    rows = []
    skipped = []
    for pid, smu, lat, lon in zip(ids, smu_ids, lats, lons):
        if smu < 0:
            skipped.append(pid)
            continue
        sub = layers[layers[cols["smu"]] == smu].copy()
        if sub.empty:
            skipped.append(pid)
            continue
        # Dominant component: highest SHARE, else lowest SEQUENCE.
        if "share" in cols:
            keep_seq = sub.loc[sub[cols["share"]].astype(float).idxmax(), cols.get("seq", cols["smu"])] \
                if "seq" in cols else None
            if "seq" in cols and keep_seq is not None:
                sub = sub[sub[cols["seq"]] == keep_seq]
        elif "seq" in cols:
            sub = sub[sub[cols["seq"]] == sub[cols["seq"]].min()]

        for _, lyr in sub.iterrows():
            bot = float(lyr[cols["bot"]])
            top = float(lyr[cols["top"]]) if "top" in cols else max(0.0, bot - 20.0)
            center = (top + bot) / 2.0
            oc = float(lyr[cols["oc"]]) if "oc" in cols and pd.notna(lyr[cols["oc"]]) else np.nan
            rows.append({
                "ID": pid, "latitude": float(lat), "longitude": float(lon),
                "depth_bottom": bot, "depth_center": center,
                "sand": float(lyr[cols["sand"]]) if pd.notna(lyr[cols["sand"]]) else np.nan,
                "clay": float(lyr[cols["clay"]]) if pd.notna(lyr[cols["clay"]]) else np.nan,
                "silt": float(lyr[cols["silt"]]) if "silt" in cols and pd.notna(lyr[cols["silt"]]) else np.nan,
                "bdod": float(lyr[cols["bulk"]]) if "bulk" in cols and pd.notna(lyr[cols["bulk"]]) else np.nan,
                "soc_pct": oc,
                "cfvo": float(lyr[cols["coarse"]]) if "coarse" in cols and pd.notna(lyr[cols["coarse"]]) else 0.0,
            })

    if not rows:
        raise RuntimeError(
            "No usable HWSD2 soil data for any point. Check the raster/DB or "
            "coordinates (all points fell on no-data / empty mapping units).")

    df = pd.DataFrame(rows)
    # Fill silt if absent: silt = 100 - sand - clay.
    miss_silt = df["silt"].isna()
    df.loc[miss_silt, "silt"] = (100.0 - df.loc[miss_silt, "sand"]
                                 - df.loc[miss_silt, "clay"]).clip(lower=0)
    df["om_pct"] = df["soc_pct"] * 1.724

    # 4. Physics + write .SOL per point (reusing the SoilGrids writer).
    physics = df.apply(
        lambda r: pd.Series(_calculate_soil_physics(
            float(r["sand"]) if pd.notna(r["sand"]) else 40.0,
            float(r["clay"]) if pd.notna(r["clay"]) else 20.0,
            float(r["om_pct"]) if pd.notna(r["om_pct"]) else 1.0,
        )), axis=1)
    df = pd.concat([df, physics], axis=1)

    log_path = os.path.join(output_sol_dir, "soil_processing_errors.log")
    with open(log_path, "w") as lf:
        lf.write(f"Log started: {datetime.now()}\n")

    success = errors = 0
    written_ids = []
    for uid in df["ID"].unique():
        subset = df[df["ID"] == uid].copy()
        try:
            _format_dssat_sol_file(subset, output_sol_dir, source_tag="HWSD2",
                                   source_name="FAO HWSD v2.0")
            success += 1
            written_ids.append(uid)
        except Exception as exc:  # noqa: BLE001
            errors += 1
            with open(log_path, "a") as lf:
                lf.write(f"ID: {uid} | Error: {exc}\n")

    if skipped:
        print(f"  [warn] {len(skipped)} point(s) had no HWSD2 mapping unit "
              f"and were skipped: {', '.join(map(str, skipped[:10]))}"
              + (" ..." if len(skipped) > 10 else ""))

    # 5. Mapping CSV (ID -> SOIL_ID == ID, matching other soil sources).
    mapping = pd.DataFrame({id_col: ids, "SOIL_ID": ids})
    mapping.to_csv(output_csv_path, index=False)

    print(f"HWSD2 processing complete. Success: {success}, Errors: {errors}, "
          f"Skipped(no-data): {len(skipped)}")
