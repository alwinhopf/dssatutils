# File: soil_soilgrids.py
# Python port of soil_soilgrids.R
#
# Reads a pre-downloaded master DSSAT .SOL file (e.g. the SoilGrids 10 km
# Harvard Dataverse files), finds the nearest soil profile for each grid
# point by great-circle distance, rewrites profile IDs, writes individual
# per-point .SOL files, and produces a mapping CSV.
#
# Reference:
#   Folberth et al. (2019) Environ. Model. Softw. 111:218-228
#   https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/1PEEY0

import os
import re
import warnings
from typing import Optional

import numpy as np
import pandas as pd

try:
    import geopandas as gpd
    from shapely.geometry import Point
    _HAS_GEOPANDAS = True
except ImportError:
    _HAS_GEOPANDAS = False
    warnings.warn("geopandas not available; falling back to haversine distance.")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _make_numeric_id(raw: str, width: int = 8) -> Optional[str]:
    """
    Strip non-digit characters from *raw* and zero-pad to *width* digits.
    Returns None if no digits remain.
    """
    digits = re.sub(r"\D+", "", str(raw))
    if not digits:
        return None
    if len(digits) >= width:
        return digits[-(width):]
    return digits.zfill(width)


def _rewrite_header_id(header_line: str, new_id: str) -> str:
    """
    Replace the first token after '*' in a profile header line with *new_id*.
    Preserves the rest of the line (spacing, lat/lon, etc.).
    """
    new_id = str(new_id)
    if len(new_id) > 10:
        new_id = new_id[-10:]
    return re.sub(r"^\*\S+", f"*{new_id}", header_line)


def _parse_lat_lon(chunk: list[str]) -> tuple[Optional[float], Optional[float]]:
    """
    Scan a small block of lines for the @SITE row and extract the first valid
    (lat, lon) pair from the following data line.
    """
    for idx, ln in enumerate(chunk):
        if ln.strip().startswith("@SITE"):
            candidates = chunk[idx + 1 : idx + 4]
            for cand in candidates:
                if not re.search(r"\d", cand):
                    continue
                nums = re.findall(r"-?\d+\.?\d*|-?\.\d+", cand)
                nums_f = [float(n) for n in nums]
                valid = [n for n in nums_f if n != -99]
                for j in range(len(valid) - 1):
                    la, lo = valid[j], valid[j + 1]
                    if abs(la) <= 90 and abs(lo) <= 180:
                        return la, lo
    return None, None


def _haversine_nearest(
    gp_lats: np.ndarray, gp_lons: np.ndarray,
    soil_lats: np.ndarray, soil_lons: np.ndarray,
) -> np.ndarray:
    """
    For each grid point return the index of the nearest soil profile using
    the haversine distance.  Vectorised over soil profiles per grid point.
    """
    R = 6_371_000.0  # Earth radius in metres
    indices = np.empty(len(gp_lats), dtype=int)
    for i, (la, lo) in enumerate(zip(gp_lats, gp_lons)):
        dlat = np.radians(soil_lats - la)
        dlon = np.radians(soil_lons - lo)
        a = (np.sin(dlat / 2) ** 2
             + np.cos(np.radians(la)) * np.cos(np.radians(soil_lats))
             * np.sin(dlon / 2) ** 2)
        indices[i] = int(np.argmin(2 * R * np.arcsin(np.sqrt(a))))
    return indices


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def process_soils_soilgrids(
    grid_points,             # GeoDataFrame or DataFrame with lon/lat columns
    source_sol_file: str,
    output_csv_path: str,
    output_sol_dir: str,
    id_col: str = "ID",
    numeric_only_ids: bool = True,
    numeric_width: int = 8,
):
    """
    Parse an external master DSSAT .SOL file, find the nearest profile for each
    grid point, write per-point .SOL files with rewritten IDs, and produce a
    mapping CSV.  Mirrors the R ``process_soils_soilgrids`` function exactly.

    Parameters
    ----------
    grid_points : GeoDataFrame | DataFrame
        Grid points with an ID column and geometry (or explicit lon/lat).
    source_sol_file : str
        Path to the master DSSAT .SOL file (one file, many profiles).
    output_csv_path : str
        Where to write the ID → SOIL_ID mapping CSV.
    output_sol_dir : str
        Directory to write individual ``<SOIL_ID>.SOL`` files.
    id_col : str
        Column in *grid_points* containing the grid-point ID.
    numeric_only_ids : bool
        If True, soil IDs are zero-padded digit strings (recommended for DSSAT).
    numeric_width : int
        Width of zero-padded numeric IDs.

    Returns
    -------
    pd.DataFrame
        Mapping table with columns: ID, SOIL_ID, SOURCE_SOIL_ID, SOIL_LAT, SOIL_LON.
    """
    print(f"Parsing external soil file: {source_sol_file}")

    if not os.path.exists(source_sol_file):
        raise FileNotFoundError(
            f"CRITICAL ERROR: The external soil file was not found at: {source_sol_file}"
        )

    os.makedirs(output_sol_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # 1. Read the master file
    # -----------------------------------------------------------------------
    with open(source_sol_file, "r", encoding="latin-1", errors="replace") as fh:
        lines = fh.readlines()

    # Remove line endings for processing; keep raw for writing
    stripped = [ln.rstrip("\n\r") for ln in lines]

    # Profile header lines start with '*'
    starts = [i for i, ln in enumerate(stripped) if ln.startswith("*")]
    if not starts:
        raise ValueError("No profiles found (no lines starting with '*').")

    ends = starts[1:] + [len(stripped)]  # exclusive end for each profile

    # Extract source IDs: first non-whitespace token after '*'
    def _extract_id(line: str) -> str:
        m = re.match(r"^\*(\S+)", line)
        return m.group(1) if m else ""

    source_ids = [_extract_id(stripped[s]) for s in starts]

    # -----------------------------------------------------------------------
    # 2. Extract coordinates from each profile block
    # -----------------------------------------------------------------------
    print("Extracting coordinates (scanning for valid Lat/Lon patterns)...")

    soil_records = []
    for i, (s, e) in enumerate(zip(starts, ends)):
        chunk_end = min(e, s + 26)
        chunk = stripped[s:chunk_end]
        lat, lon = _parse_lat_lon(chunk)
        if lat is None or lon is None:
            continue
        soil_records.append({
            "source_soil_id": source_ids[i],
            "lat": lat,
            "lon": lon,
            "start_line": s,
            "end_line": e - 1,
        })

    soil_df = pd.DataFrame(soil_records).dropna(subset=["lat", "lon"])
    if soil_df.empty:
        raise ValueError(
            "No soil profiles with valid lat/lon were found in the master .SOL. "
            "Check the @SITE blocks."
        )

    # -----------------------------------------------------------------------
    # 3. Prepare grid-point coordinates
    # -----------------------------------------------------------------------
    if _HAS_GEOPANDAS and hasattr(grid_points, "geometry"):
        gdf = grid_points.to_crs("EPSG:4326")
        gp_lons = gdf.geometry.x.values
        gp_lats = gdf.geometry.y.values
        gp_df = gdf
    else:
        gp_df = pd.DataFrame(grid_points)
        # Try common coordinate column names
        lon_candidates = ["lon", "Lon", "LON", "x", "X", "LONG", "longitude", "Longitude"]
        lat_candidates = ["lat", "Lat", "LAT", "y", "Y", "LATI", "latitude", "Latitude"]
        lon_col = next((c for c in lon_candidates if c in gp_df.columns), None)
        lat_col = next((c for c in lat_candidates if c in gp_df.columns), None)
        if lon_col is None or lat_col is None:
            raise ValueError(
                "grid_points must be a GeoDataFrame, or a DataFrame with lon/lat columns."
            )
        gp_lons = gp_df[lon_col].values.astype(float)
        gp_lats = gp_df[lat_col].values.astype(float)

    if id_col not in gp_df.columns:
        raise ValueError(f"id_col '{id_col}' not found in grid_points.")

    # -----------------------------------------------------------------------
    # 4. Nearest soil profile for each grid point
    # -----------------------------------------------------------------------
    soil_lats = soil_df["lat"].values
    soil_lons = soil_df["lon"].values

    if _HAS_GEOPANDAS:
        soil_gdf = gpd.GeoDataFrame(
            soil_df,
            geometry=gpd.points_from_xy(soil_lons, soil_lats),
            crs="EPSG:4326",
        )
        gp_gdf = gpd.GeoDataFrame(
            gp_df.reset_index(drop=True),
            geometry=gpd.points_from_xy(gp_lons, gp_lats),
            crs="EPSG:4326",
        )
        nearest_idx = gp_gdf.geometry.apply(
            lambda pt: soil_gdf.geometry.distance(pt).idxmin()
        ).values
    else:
        nearest_idx = _haversine_nearest(gp_lats, gp_lons, soil_lats, soil_lons)

    nearest = soil_df.iloc[nearest_idx].reset_index(drop=True)

    point_ids_raw = gp_df[id_col].astype(str).tolist()
    if numeric_only_ids:
        soil_ids = [_make_numeric_id(pid, numeric_width) for pid in point_ids_raw]
        # Fallback for any None
        for i, sid in enumerate(soil_ids):
            if sid is None:
                soil_ids[i] = str(i + 1).zfill(numeric_width)
    else:
        soil_ids = point_ids_raw

    mapping = pd.DataFrame({
        "ID": point_ids_raw,
        "SOIL_ID": soil_ids,
        "SOURCE_SOIL_ID": nearest["source_soil_id"].values,
        "SOIL_LAT": nearest["lat"].values,
        "SOIL_LON": nearest["lon"].values,
        "START_LINE": nearest["start_line"].values,
        "END_LINE": nearest["end_line"].values,
    })

    # -----------------------------------------------------------------------
    # 5. Write mapping CSV
    # -----------------------------------------------------------------------
    mapping[["ID", "SOIL_ID", "SOURCE_SOIL_ID", "SOIL_LAT", "SOIL_LON"]].to_csv(
        output_csv_path, index=False
    )

    # -----------------------------------------------------------------------
    # 6. Write per-point .SOL files with rewritten profile IDs
    # -----------------------------------------------------------------------
    print("Writing individual .SOL files...")
    for _, mrow in mapping.iterrows():
        s = int(mrow["START_LINE"])
        e = int(mrow["END_LINE"]) + 1  # inclusive → exclusive
        prof_lines = stripped[s:e]

        # Rewrite the header ID
        prof_lines[0] = _rewrite_header_id(prof_lines[0], mrow["SOIL_ID"])

        out_file = os.path.join(output_sol_dir, f"{mrow['SOIL_ID']}.SOL")
        with open(out_file, "w", encoding="latin-1") as fh:
            fh.write("\n".join(prof_lines) + "\n")

    result = mapping[["ID", "SOIL_ID", "SOURCE_SOIL_ID", "SOIL_LAT", "SOIL_LON"]]
    print(f"Done. {len(mapping)} soil profiles written to: {output_sol_dir}")
    return result
