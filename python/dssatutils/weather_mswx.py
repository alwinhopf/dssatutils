# Weather source: MSWX global meteorological forcing (local NetCDF cache).

from .weather_gridded_common import process_local_netcdf_weather

_VARS = {
    "TMAX": {"tokens": ["tmax", "tasmax"], "kind": "temp", "required": True},
    "TMIN": {"tokens": ["tmin", "tasmin"], "kind": "temp", "required": True},
    "TMEAN": {"tokens": ["temp", "tas"], "kind": "temp"},
    "RAIN": {"tokens": ["precip", "pr", "mswx"], "kind": "rain", "required": True},
    "SRAD": {"tokens": ["srad", "rsds", "shortwave"], "kind": "srad", "required": True},
    "WIND": {"tokens": ["wind", "sfcwind"], "kind": "wind"},
    "RH2M": {"tokens": ["rh", "hurs"], "kind": "rh"},
}


def process_weather_mswx(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    mswx_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local MSWX NetCDF files."""
    print(f"--- Starting MSWX Processing (Years: {start_year}-{end_year}) ---")
    written = process_local_netcdf_weather(
        shapefile, int(start_year), int(end_year), output_dir,
        id_col, lat_col, lon_col, log_file, mswx_nc_dir, _VARS,
        "MSWX", "MSWX")
    print(f"\nMSWX processing complete: {written} point(s) written.\n")
