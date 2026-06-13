# dssatutils

Shared **weather** and **soil** download utilities for DSSAT gridded / spatial
crop-model pipelines. One versioned home for the download logic that was
previously duplicated between:

- **DSSAT_Gridded_Run_Tutorial** (source of truth — R + Python)
- **DSSAT_ML_Phenology_Prediction** (R)

Each function fetches data from a public source and writes DSSAT-format
`.WTH` (weather) or `.SOL` (soil) files for a set of grid points.

> **Private repository.** Install requires access to `github.com/alwinhopf/dssatutils`.

## What's inside

| Domain | Sources (function name is the same in R and Python) |
|---|---|
| Weather | `process_weather_daymet`, `process_weather_gridmet`, `process_weather_nasapower`, `process_weather_openmeteo`, `process_weather_agera5`, `process_weather_nasapower_chirps` |
| Soil | `process_soils_ssurgo`, `process_soils_ssurgo_alderman`, `process_soils_soilgrids`, `process_soils_soilgrids_online`, `process_soils_hwsd` |

Coverage notes: Daymet = North America; GridMET/SSURGO = USA; NASA POWER /
Open-Meteo / AgERA5 / SoilGrids / HWSD2 = global. AgERA5 needs a free
Copernicus CDS API key; CHIRPS fuses NASA POWER with high-res rainfall (50S–50N).

## Install

### R
```r
# install.packages("remotes")
remotes::install_github("alwinhopf/dssatutils@v0.1.0")
library(dssatutils)
```

### Python
```bash
pip install "git+https://github.com/alwinhopf/dssatutils.git@v0.1.0"
# AgERA5 backend (optional, needs a Copernicus CDS key):
pip install "dssatutils[agera5] @ git+https://github.com/alwinhopf/dssatutils.git@v0.1.0"
```
or pin in `requirements.txt`:
```
dssatutils @ git+https://github.com/alwinhopf/dssatutils.git@v0.1.0
```

```python
from dssatutils import process_weather_nasapower, process_soils_ssurgo
```

## Versioning

Semantic versioning with Git tags. **Consumer repos always pin to a tag**
(`@v0.1.0`), never `main`, so upstream changes never break a pipeline until you
deliberately bump the pin. Workflow: branch → CI smoke tests → merge → tag
`vX.Y.Z` → bump the pin in each consumer repo.

## Known limitations / notes

- **GridMET** RH2M and TDEW are *estimated* (`TDEW ≈ TMIN − 2.5`, RH from the
  diurnal temperature range), not measured.
- **Open-Meteo** has no daily dewpoint/RH, so it writes `TDEW = -99`, `RH2M = -99`
  (DSSAT-valid missing). ET methods needing RH will degrade.
- **TAV/AMP** is computed via `DSSAT::calc_TAV/calc_AMP` for GridMET but hand-rolled
  (monthly-mean amplitude) for the other sources — values are close but not
  identical. Consolidating into one shared helper is a planned cleanup.
- **SoilGrids online** mode is selected via a package-level `USE_REST_API` flag
  (REST vs VRT); a future version will make it a function argument.

See `SHARED_UTILS_MIGRATION.md` in the Gridded Run Tutorial repo for the full
extraction history and the remaining packaging-polish checklist.
