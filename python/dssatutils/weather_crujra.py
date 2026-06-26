# Weather source: CRU-JRA bias-corrected reanalysis (local NetCDF cache).

from .weather_gridded_common import process_local_netcdf_weather

_VARS = {
    "TMAX": {"tokens": ["tmax", "tasmax"], "kind": "temp", "required": True},
    "TMIN": {"tokens": ["tmin", "tasmin"], "kind": "temp", "required": True},
    "TMEAN": {"tokens": ["tas", "tmean"], "kind": "temp"},
    "RAIN": {"tokens": ["pr", "precip"], "kind": "rain", "required": True},
    "SRAD": {"tokens": ["rsds", "srad"], "kind": "srad", "required": True},
    "WIND": {"tokens": ["wind", "sfcwind"], "kind": "wind"},
    "RH2M": {"tokens": ["rh", "hurs"], "kind": "rh"},
}


def process_weather_crujra(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    crujra_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local CRU-JRA NetCDF files."""
    print(f"--- Starting CRU-JRA Processing (Years: {start_year}-{end_year}) ---")
    written = process_local_netcdf_weather(
        shapefile, int(start_year), int(end_year), output_dir,
        id_col, lat_col, lon_col, log_file, crujra_nc_dir, _VARS,
        "CRU-JRA", "CRUJ")
    print(f"\nCRU-JRA processing complete: {written} point(s) written.\n")
