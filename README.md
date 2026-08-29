# dssatutils

> **AI agents & maintainers:** read [`../AGENTS.md`](../AGENTS.md) before editing this repo.

Shared **weather** and **soil** download utilities for DSSAT gridded / spatial
crop-model pipelines. One versioned home for the download logic that was
previously duplicated between:

- **DSSAT_Gridded_Run_Tutorial** (source of truth — R + Python)
- **DSSAT_ML_Phenology_Prediction** (R)

Each function fetches data from a public source and writes DSSAT-format
`.WTH` (weather) or `.SOL` (soil) files for a set of grid points.

> **GitHub install.** Install requires access to `github.com/alwinhopf/dssatutils`.
> If Git prompts for authentication, configure Git Credential Manager, SSH keys,
> or a GitHub token before running the install command.

## What's inside

| Domain | Sources (function name is the same in R and Python) |
|---|---|
| Weather | `process_weather_daymet`, `process_weather_gridmet`, `process_weather_nasapower`, `process_weather_openmeteo`, `process_weather_agera5`, `process_weather_nasapower_chirps`, `process_weather_nasapower_chirps_v3`, `extract_chirps_v3_rainfall`, `merge_rainfall_into_weather`, `process_weather_cmfd`, `process_weather_dwd`, `process_weather_eobs`, `process_weather_xavier`, `process_weather_era5_land`, `process_weather_chelsa_w5e5`, `process_weather_agmerra`, `process_weather_agcfsr`, `process_weather_silo`, `process_weather_prism`, `process_weather_mswx`, `process_weather_mswep`, `process_weather_crujra`, `process_weather_terraclimate`, `process_weather_aphrodite`, `process_weather_anusplin`, `process_weather_tamsat`, `process_weather_ghcn`, `process_weather_pgf`, `process_weather_merra2` |
| Soil | `process_soils_ssurgo`, `process_soils_ssurgo_alderman`, `process_soils_polaris`, `process_soils_soilgrids`, `process_soils_soilgrids_online`, `process_soils_hwsd`, `process_soils_agmip`, `process_soils_hihydrosoil`, `process_soils_slga`, `process_soils_wise30sec`, `process_soils_wosis`, `process_soils_gnatsgo`, `process_soils_isdasoil`, `process_soils_lucas`, `process_soils_gsde`, `process_soils_china`, `process_soils_febr`, `process_soils_slc`, `process_soils_esdb`, `process_soils_openlandmap` |

Coverage notes: Daymet = North America; GridMET/SSURGO/gNATSGO/POLARIS = USA; NASA POWER /
Open-Meteo / AgERA5 / ERA5-Land / SoilGrids / HWSD2 = global; iSDAsoil = Africa;
LUCAS = Europe topsoil; CMFD = China; DWD = Germany; E-OBS = Europe; Xavier = Brazil;
SILO/SLGA = Australia; PRISM = CONUS; CHELSA-W5E5 / AgMERRA / AgCFSR / HiHydroSoil /
MSWX / MSWEP / CRU-JRA / TerraClimate / WISE30sec / WoSIS = global or near-global.
AgMIP/Han = global 5 arc-min DSSAT-ready country `.SOL` files (local download required).
Newer regional fills: APHRODITE = monsoon Asia rainfall (NASA-POWER hybrid); ANUSPLIN = Canada temperature/precipitation only. ANUSPLIN is intentionally rejected as standalone DSSAT forcing unless an SRAD layer is also supplied, because its core product cannot provide a physically complete WTH file;
TAMSAT = Africa rainfall (NASA-POWER hybrid); PGF / MERRA-2 = global reanalysis; GHCN-Daily =
global station obs (live NOAA download, nearest-station). Soil: GSDE = global 1 km 8-layer;
China BNU = China; FEBR/Embrapa = Brazil; SLC = Canada; ESDB = Europe full profile;
OpenLandMap = global 250 m (live COG sampling, no local data).
AgERA5, ERA5-Land, and E-OBS require a free Copernicus CDS API key; CHIRPS v2 fuses
NASA POWER with high-res rainfall (50S-50N). CHIRPS v3 is available via
`extract_chirps_v3_rainfall()` and `process_weather_nasapower_chirps_v3()` with
`rnl` (ERA5 daily disaggregation, full historical period) and `sat` (IMERG daily
disaggregation, recent period) options; v3 coverage is 60S-60N and currently uses
p05 daily NetCDFs. The v3 helper defaults to monthly NetCDF caching because
yearly v3 daily files are roughly 23.5 GiB/product-year, while a typical monthly
file is a few hundred MiB. POLARIS is a 30 m probabilistic
disaggregation of SSURGO (Chaney et al. 2019); `process_soils_polaris` builds
water limits from its van Genuchten curve and takes a `stat` argument (default
`"p50"`, the deterministic median — p5/p95 percentile layers are reserved for a
future uncertainty-ensemble layer).

Live CHIRPS validation is available but skipped by default because it downloads
real NetCDFs. Run it explicitly with
`DSSATUTILS_RUN_LIVE_CHIRPS=1 python -m pytest tests/test_chirps_live.py -m live -q -s`.
The test caches data in `.live_cache/`, generates NASA POWER + CHIRPS v2 and
NASA POWER + CHIRPS v3 DSSAT `.WTH` files for a real point, and compares daily
rainfall over the overlapping period.

## Install

### R
```r
# install.packages("remotes")
remotes::install_github("alwinhopf/dssatutils@e9c859fa1d915623df23e2eb13084cb085dbfe3e")
library(dssatutils)
```

### Python
```bash
pip install "git+https://github.com/alwinhopf/dssatutils.git@e9c859fa1d915623df23e2eb13084cb085dbfe3e"
# CDS-backed weather sources: AgERA5, ERA5-Land, optional E-OBS CDS mode.
pip install "dssatutils[cds] @ git+https://github.com/alwinhopf/dssatutils.git@e9c859fa1d915623df23e2eb13084cb085dbfe3e"
```
or pin in `requirements.txt`:
```
dssatutils @ git+https://github.com/alwinhopf/dssatutils.git@e9c859fa1d915623df23e2eb13084cb085dbfe3e
```

```python
from dssatutils import process_weather_nasapower, process_soils_ssurgo
```

## Credentials

Most sources are keyless or use local files you download separately. Copernicus
CDS-backed sources (`process_weather_agera5`, `process_weather_era5_land`, and
`process_weather_eobs(..., eobs_use_cds=TRUE/True)`) need a free CDS Personal
Access Token and accepted dataset licences.

R:
```r
library(dssatutils)
setup_cds_credentials()
```

Python:
```python
from dssatutils import setup_cds_credentials
setup_cds_credentials()
```

The helper uses `CDSAPI_KEY`/`CDSAPI_URL`, imports an existing `~/.cdsapirc`, or
prompts in an interactive session. It writes a cdsapi-compatible `.cdsapirc`;
the R helper also stores the token for `ecmwfr`.

## Configuration

Package-level defaults live in `config.yml` and are read by both R and Python.
The installed package carries the same defaults in `inst/config.yml` and
`python/dssatutils/config.yml`. Callers can override those defaults by setting
`DSSATUTILS_CONFIG` to another YAML file or by passing explicit function
arguments. Consumer pipeline `config.yml` files remain the study-level source of
truth and are merged over the package defaults.

Currently configured package defaults include CDS URL, Open-Meteo rate-limit
settings, CHIRPS v3 product/cache/download settings, and SoilGrids Online
REST/VRT behavior.

For AgERA5, consumer pipelines should prefer `agera5_backend = "timeseries"`.
The time-series backend requests all seven DSSAT weather variables together and
stores CSVs by year and globally anchored AgERA5 grid chunk. With
`agera5_timeseries_chunk_degrees = 0.1`, each cache entry represents one
canonical 0.1-degree AgERA5 cell, so crops, soils, point subsets, and model-grid
resolutions that select the same cell reuse the same download. Larger values
group cells into fixed global tiles and reduce request count at the cost of more
downloaded data. The legacy `gridded` backend remains available for callers that
need the original daily-NetCDF ZIPs.

## Versioning

Semantic versioning with Git tags. **Consumer repos always pin to a tag**
(`@vX.Y.Z`), never `main`, so upstream changes never break a pipeline until you
deliberately bump the pin. Workflow: branch → CI smoke tests → merge → tag
`vX.Y.Z` → bump the pin in each consumer repo.

## Known limitations / notes

- **Optional weather repair and QA** is available after provider downloads via
  `repair_weather_missing_values()`, `repair_weather_date_gaps()`,
  `repair_weather_temperature_inversions()`, and `audit_weather_quality()`.
  Repair functions only modify short runs with valid before/after neighbors;
  the audit writes flag-only findings to CSV and appends notes to
  `weather_repair.log`. Provider `NA`/`NaN`/infinite values are normalized to
  DSSAT's numeric `-99` marker before writing, and the repair reader also accepts
  those literal tokens in older cached files.
- **Weather-file validation** via `is_wth_valid()` understands DSSAT's
  fixed-width daily rows (including adjacent negative fields), requires
  consecutive dates, and rejects physically impossible forcing while retaining
  the standard `-99` missing-value sentinel. Callers can pass
  `required_columns` to reject `-99` in model-essential fields while still
  permitting optional missing humidity/wind inputs. AgERA5 applies the same
  physical checks before writing. Its recommended time-series cache uses globally
  anchored 0.1-degree cells or fixed tiles; the legacy gridded cache remains
  keyed by the exact requested geographic area.
- **GridMET** RH2M and TDEW are *estimated* (`TDEW ≈ TMIN − 2.5`, RH from the
  diurnal temperature range), not measured.
- **Open-Meteo** uses the API's `dew_point_2m_mean` and
  `relative_humidity_2m_mean` for DSSAT `TDEW` and `RH2M`. The default
  `era5_seamless` model combines ERA5-Land temperature/humidity with ERA5
  forcing fields so radiation, precipitation, and wind remain complete. The R
  adapter runs its rate-limited request stream sequentially in-process.
- **TAV/AMP** is computed via `DSSAT::calc_TAV/calc_AMP` for GridMET but hand-rolled
  (monthly-mean amplitude) for the other sources — values are close but not
  identical. Consolidating into one shared helper is a planned cleanup.
- **SoilGrids Online** defaults are controlled in `config.yml`
  (`soil.soilgrids_online.use_rest_api`; `false` = VRT, `true` = REST). The
  gridded tutorial exposes this through its `soilgrids_mode` key.

See `SHARED_UTILS_MIGRATION.md` in the Gridded Run Tutorial repo for the full
extraction history and the remaining packaging-polish checklist.
