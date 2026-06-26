# Weather sources: AgMERRA and AgCFSR AgMIP climate forcing datasets.

from .weather_gridded_common import process_local_netcdf_weather

_VARS = {
    "TMAX": {"tokens": ["tmax", "tasmax"], "aliases": ["tmax", "tasmax"], "kind": "temp", "required": True},
    "TMIN": {"tokens": ["tmin", "tasmin"], "aliases": ["tmin", "tasmin"], "kind": "temp", "required": True},
    "TMEAN": {"tokens": ["tmean", "tas"], "aliases": ["tmean", "tas"], "kind": "temp"},
    "RAIN": {"tokens": ["prate", "precip", "pr"], "aliases": ["prate", "precip", "pr"], "kind": "rain", "required": True},
    "SRAD": {"tokens": ["srad", "rsds", "swdown"], "aliases": ["srad", "rsds", "swdown"], "kind": "srad", "required": True},
    "WIND": {"tokens": ["wind", "sfcwind", "wnd"], "aliases": ["wind", "sfcwind", "wnd"], "kind": "wind"},
    "TDEW": {"tokens": ["tdew", "dew"], "aliases": ["tdew", "dewpoint"], "kind": "temp"},
    "RH2M": {"tokens": ["rh", "rhum", "hurs"], "aliases": ["rh", "rhum", "hurs"], "kind": "rh"},
}


def _process(shapefile, start_year, end_year, output_dir, id_col, lat_col,
             lon_col, log_file, nc_dir, product, insi):
    start_year = max(int(start_year), 1980)
    end_year = min(int(end_year), 2010)
    print(f"--- Starting {product} Processing (Years: {start_year}-{end_year}) ---")
    written = process_local_netcdf_weather(
        shapefile, start_year, end_year, output_dir,
        id_col, lat_col, lon_col, log_file, nc_dir, _VARS,
        product, insi, refht=2.0, wndht=2.0)
    print(f"\n{product} processing complete: {written} point(s) written.\n")


def process_weather_agmerra(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    agmerra_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local AgMERRA NetCDF files (1980-2010)."""
    _process(shapefile, start_year, end_year, output_dir, id_col, lat_col,
             lon_col, log_file, agmerra_nc_dir, "AgMERRA", "AGMR")


def process_weather_agcfsr(
    shapefile, start_year, end_year, output_dir,
    id_col, lat_col, lon_col, n_cores, log_file,
    agcfsr_nc_dir: str,
) -> None:
    """Write DSSAT .WTH files from local AgCFSR NetCDF files (1980-2010)."""
    _process(shapefile, start_year, end_year, output_dir, id_col, lat_col,
             lon_col, log_file, agcfsr_nc_dir, "AgCFSR", "AGCF")
