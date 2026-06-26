# Weather source: CHELSA-W5E5 daily (global, 30 arcsec, 1979-2016).

from .weather_gridded_common import process_local_netcdf_weather

_VARS = {
    "TMAX": {"tokens": ["tasmax"], "aliases": ["tasmax"], "kind": "temp", "required": True},
    "TMIN": {"tokens": ["tasmin"], "aliases": ["tasmin"], "kind": "temp", "required": True},
    "TMEAN": {"tokens": ["tas_", "tas."], "aliases": ["tas"], "kind": "temp"},
    "RAIN": {"tokens": ["pr"], "aliases": ["pr", "precip"], "kind": "rain", "required": True},
    "SRAD": {"tokens": ["rsds"], "aliases": ["rsds"], "kind": "srad", "required": True},
}


def process_weather_chelsa_w5e5(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    chelsa_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local CHELSA-W5E5 daily NetCDF files.

    Expected files contain tasmax, tasmin, pr, and rsds (one variable per file is
    fine). The dataset ends in 2016; requested later years are clipped.
    """
    end_year = min(int(end_year), 2016)
    print(f"--- Starting CHELSA-W5E5 Processing (Years: {start_year}-{end_year}) ---")
    written = process_local_netcdf_weather(
        shapefile, int(start_year), end_year, output_dir,
        id_col, lat_col, lon_col, log_file, chelsa_nc_dir, _VARS,
        "CHELSA-W5E5", "CHW5", refht=2.0, wndht=2.0)
    print(f"\nCHELSA-W5E5 processing complete: {written} point(s) written.\n")
