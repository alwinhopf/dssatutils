# Weather source: NASA MERRA-2 reanalysis daily aggregates (local NetCDF).
#
# Global ~0.5 x 0.625 degree reanalysis (1980-present). Read locally-prepared
# daily NetCDF (e.g. statM/aggregated T2MMAX/T2MMIN/PRECTOT/SWGDN/SPEED) through
# the shared extractor. Temperatures arrive in Kelvin and are auto-converted.

from .weather_gridded_common import process_local_netcdf_weather

_VARS = {
    "TMAX": {"tokens": ["t2mmax", "tmax", "tasmax"], "kind": "temp", "required": True},
    "TMIN": {"tokens": ["t2mmin", "tmin", "tasmin"], "kind": "temp", "required": True},
    "TMEAN": {"tokens": ["t2mmean", "t2m", "tas"], "kind": "temp"},
    "RAIN": {"tokens": ["prectot", "precip", "pr"], "kind": "rain", "required": True},
    "SRAD": {"tokens": ["swgdn", "swgnt", "rsds", "srad"], "kind": "srad", "required": True},
    "WIND": {"tokens": ["speed", "wind", "sfcwind"], "kind": "wind"},
    "RH2M": {"tokens": ["rh2m", "rh", "hurs"], "kind": "rh"},
}


def process_weather_merra2(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    merra2_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local MERRA-2 daily NetCDF files (1980-present)."""
    print(f"--- Starting MERRA-2 Processing (Years: {start_year}-{end_year}) ---")
    written = process_local_netcdf_weather(
        shapefile, max(int(start_year), 1980), int(end_year), output_dir,
        id_col, lat_col, lon_col, log_file, merra2_nc_dir, _VARS,
        "NASA MERRA-2", "MER2")
    print(f"\nMERRA-2 processing complete: {written} point(s) written.\n")
