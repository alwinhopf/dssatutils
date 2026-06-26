# Weather source: Princeton Global Forcing (PGF) daily reanalysis (local NetCDF).
#
# Global 0.25 degree daily forcing (1948-present in the v3 line). A full-variable
# alternative/gap-fill to AgMERRA/CRU-JRA via the shared NetCDF extractor.

from .weather_gridded_common import process_local_netcdf_weather

_VARS = {
    "TMAX": {"tokens": ["tmax", "tasmax"], "kind": "temp", "required": True},
    "TMIN": {"tokens": ["tmin", "tasmin"], "kind": "temp", "required": True},
    "TMEAN": {"tokens": ["tas", "tmean"], "kind": "temp"},
    "RAIN": {"tokens": ["prcp", "precip", "pr"], "kind": "rain", "required": True},
    "SRAD": {"tokens": ["dswrf", "rsds", "srad", "swdown"], "kind": "srad", "required": True},
    "WIND": {"tokens": ["wind", "wnd", "sfcwind"], "kind": "wind"},
    "RH2M": {"tokens": ["rh", "hurs"], "kind": "rh"},
}


def process_weather_pgf(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    pgf_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local Princeton Global Forcing NetCDF files."""
    print(f"--- Starting PGF Processing (Years: {start_year}-{end_year}) ---")
    written = process_local_netcdf_weather(
        shapefile, int(start_year), int(end_year), output_dir,
        id_col, lat_col, lon_col, log_file, pgf_nc_dir, _VARS,
        "Princeton Global Forcing", "PGF")
    print(f"\nPGF processing complete: {written} point(s) written.\n")
