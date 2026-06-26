# Weather source: SILO Australian daily gridded climate data.

from .weather_gridded_common import process_local_netcdf_weather

_VARS = {
    "TMAX": {"tokens": ["max_temp", "tmax", "maxt"], "aliases": ["max_temp", "tmax", "maxt"], "kind": "temp", "required": True},
    "TMIN": {"tokens": ["min_temp", "tmin", "mint"], "aliases": ["min_temp", "tmin", "mint"], "kind": "temp", "required": True},
    "RAIN": {"tokens": ["rain", "precip", "ppt"], "aliases": ["rain", "precip", "ppt"], "kind": "rain", "required": True},
    "SRAD": {"tokens": ["radiation", "srad", "solar"], "aliases": ["radiation", "srad", "solar"], "kind": "srad", "required": True},
    # SILO distributes vapour pressure ("vp", hPa), not dewpoint: convert via kind="vp".
    "TDEW": {"tokens": ["vp"], "aliases": ["vp", "vapour", "vapor"], "kind": "vp"},
    "WIND": {"tokens": ["wind"], "aliases": ["wind", "sfcwind"], "kind": "wind"},
}


def process_weather_silo(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    silo_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local SILO NetCDF files for Australia."""
    print(f"--- Starting SILO Processing (Years: {start_year}-{end_year}) ---")
    written = process_local_netcdf_weather(
        shapefile, int(start_year), int(end_year), output_dir,
        id_col, lat_col, lon_col, log_file, silo_nc_dir, _VARS,
        "SILO Australia", "SILO", refht=2.0, wndht=2.0)
    print(f"\nSILO processing complete: {written} point(s) written.\n")
