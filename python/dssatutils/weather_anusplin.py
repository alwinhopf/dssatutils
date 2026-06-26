# Weather source: ANUSPLIN (Agriculture & Agri-Food Canada) daily 10 km grids.
#
# Station-spline interpolation covering ALL of Canada (1950-2015), filling the
# gap above Daymet's ~52N ceiling. Distributed as local NetCDF (maxt/mint/pcp);
# SRAD/RH/wind are not in the core product and are written DSSAT-missing (-99).

from .weather_gridded_common import process_local_netcdf_weather

_VARS = {
    "TMAX": {"tokens": ["maxt", "tmax", "tasmax"], "kind": "temp", "required": True},
    "TMIN": {"tokens": ["mint", "tmin", "tasmin"], "kind": "temp", "required": True},
    "RAIN": {"tokens": ["pcp", "precip", "pr", "rain"], "kind": "rain", "required": True},
}


def process_weather_anusplin(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    anusplin_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local ANUSPLIN (Canada) NetCDF files."""
    end_year = min(int(end_year), 2015)
    print(f"--- Starting ANUSPLIN Processing (Years: {start_year}-{end_year}) ---")
    written = process_local_netcdf_weather(
        shapefile, int(start_year), end_year, output_dir,
        id_col, lat_col, lon_col, log_file, anusplin_nc_dir, _VARS,
        "ANUSPLIN Canada", "ANUS")
    print(f"\nANUSPLIN processing complete: {written} point(s) written.\n")
