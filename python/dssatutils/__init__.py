"""dssatutils — shared weather & soil download utilities for DSSAT pipelines.

Public API: the ``process_*`` functions below write DSSAT ``.WTH`` (weather) and
``.SOL`` (soil) files for a set of grid points. Underscore-prefixed helpers in the
submodules are internal and not re-exported here.

Each submodule pulls in heavy geospatial dependencies (xarray, rasterio,
geopandas, ...) only when imported. To keep ``import dssatutils`` cheap and to
avoid hard-failing when an optional backend (e.g. ``cdsapi`` for AgERA5) is
missing, submodules are imported lazily via ``__getattr__`` (PEP 562).
"""
from importlib import import_module

__version__ = "0.4.0"

# public name -> submodule that defines it
_EXPORTS = {
    "process_weather_daymet": "weather_daymet",
    "process_weather_gridmet": "weather_gridmet",
    "process_weather_nasapower": "weather_nasapower",
    "process_weather_agera5": "weather_agera5",
    "process_weather_openmeteo": "weather_openmeteo",
    "process_weather_nasapower_chirps": "weather_nasapower_chirps",
    "extract_chirps_v3_rainfall": "weather_chirps_v3",
    "process_weather_nasapower_chirps_v3": "weather_chirps_v3",
    "merge_rainfall_into_weather": "weather_rainfall_merge",
    "process_weather_cmfd": "weather_cmfd",
    "process_weather_dwd": "weather_dwd",
    "process_weather_eobs": "weather_eobs",
    "process_weather_xavier": "weather_xavier",
    "process_weather_era5_land": "weather_era5land",
    "process_weather_chelsa_w5e5": "weather_chelsa_w5e5",
    "process_weather_agmerra": "weather_agmip",
    "process_weather_agcfsr": "weather_agmip",
    "process_weather_silo": "weather_silo",
    "process_weather_prism": "weather_prism",
    "process_weather_mswx": "weather_mswx",
    "process_weather_mswep": "weather_mswep",
    "process_weather_crujra": "weather_crujra",
    "process_weather_terraclimate": "weather_terraclimate",
    "process_weather_aphrodite": "weather_aphrodite",
    "process_weather_anusplin": "weather_anusplin",
    "process_weather_tamsat": "weather_tamsat",
    "process_weather_ghcn": "weather_ghcn",
    "process_weather_pgf": "weather_pgf",
    "process_weather_merra2": "weather_merra2",
    "setup_cds_credentials": "credentials",
    "era5land_set_cds_key": "credentials",
    "repair_weather_missing_values": "weather_repair",
    "repair_weather_file_missing_values": "weather_repair",
    "repair_weather_temperature_inversions": "weather_repair",
    "repair_weather_file_temperature_inversions": "weather_repair",
    "repair_weather_date_gaps": "weather_repair",
    "repair_weather_file_date_gaps": "weather_repair",
    "audit_weather_quality": "weather_repair",
    "audit_weather_file_quality": "weather_repair",
    "process_soils_soilgrids": "soil_soilgrids",
    "process_soils_soilgrids_online": "soil_soilgrids_online",
    "process_soils_ssurgo": "soil_ssurgo",
    "process_soils_ssurgo_alderman": "soil_ssurgo_alderman",
    "pull_profile_by_name_alderman": "soil_ssurgo_alderman",
    "pull_profile_by_coords_alderman": "soil_ssurgo_alderman",
    "process_soils_polaris": "soil_polaris",
    "process_soils_hwsd": "soil_hwsd",
    "process_soils_agmip": "soil_agmip",
    "process_soils_hihydrosoil": "soil_hihydrosoil",
    "process_soils_slga": "soil_slga",
    "process_soils_wise30sec": "soil_wise30sec",
    "process_soils_wosis": "soil_wosis",
    "process_soils_gsde": "soil_gsde",
    "process_soils_china": "soil_china",
    "process_soils_febr": "soil_febr",
    "process_soils_slc": "soil_slc",
    "process_soils_esdb": "soil_esdb",
    "process_soils_openlandmap": "soil_openlandmap",
    "process_soils_gnatsgo": "soil_gnatsgo",
    "process_soils_isdasoil": "soil_isdasoil",
    "process_soils_lucas": "soil_lucas",
}

__all__ = list(_EXPORTS)


def __getattr__(name):  # PEP 562 lazy attribute loading
    mod_name = _EXPORTS.get(name)
    if mod_name is not None:
        module = import_module(f".{mod_name}", __name__)
        return getattr(module, name)
    if name.startswith(("weather_", "soil_")):
        try:
            module = import_module(f".{name}", __name__)
        except ModuleNotFoundError as exc:
            if exc.name == f"{__name__}.{name}":
                raise AttributeError(f"module 'dssatutils' has no attribute {name!r}") from exc
            raise
        globals()[name] = module
        return module
    raise AttributeError(f"module 'dssatutils' has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + __all__)
