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

__version__ = "0.1.0"

# public name -> submodule that defines it
_EXPORTS = {
    "process_weather_daymet": "weather_daymet",
    "process_weather_gridmet": "weather_gridmet",
    "process_weather_nasapower": "weather_nasapower",
    "process_weather_agera5": "weather_agera5",
    "process_weather_openmeteo": "weather_openmeteo",
    "process_weather_nasapower_chirps": "weather_nasapower_chirps",
    "process_soils_soilgrids": "soil_soilgrids",
    "process_soils_soilgrids_online": "soil_soilgrids_online",
    "process_soils_ssurgo": "soil_ssurgo",
    "process_soils_ssurgo_alderman": "soil_ssurgo_alderman",
    "process_soils_hwsd": "soil_hwsd",
}

__all__ = list(_EXPORTS)


def __getattr__(name):  # PEP 562 lazy attribute loading
    mod_name = _EXPORTS.get(name)
    if mod_name is None:
        raise AttributeError(f"module 'dssatutils' has no attribute {name!r}")
    module = import_module(f".{mod_name}", __name__)
    return getattr(module, name)


def __dir__():
    return sorted(list(globals()) + __all__)
