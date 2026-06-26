# Soil source: FEBR / Embrapa Brazilian soil profiles (processed CSV).
#
# FEBR ("Free Brazilian Repository for Open Soil Data") and Embrapa profiles give
# Brazil-specific soil that pairs with the BR-DWGD/Xavier weather backend. This
# backend ingests a harmonized per-layer CSV and selects the nearest profile per
# grid point (same mechanism as WoSIS). Required CSV columns: profile_id,
# latitude, longitude, depth_bottom, sand, clay, silt, bdod, soc_pct.

from .soil_wosis import process_soils_wosis


def process_soils_febr(
    grid_points, febr_profile_csv: str, output_csv_path: str, output_sol_dir: str,
    id_col: str = "ID", lat_col: str = "LAT", long_col: str = "LONG",
) -> None:
    """Build DSSAT .SOL files from a harmonized FEBR/Embrapa Brazil layer CSV."""
    return process_soils_wosis(
        grid_points, febr_profile_csv, output_csv_path, output_sol_dir,
        id_col=id_col, lat_col=lat_col, long_col=long_col)
