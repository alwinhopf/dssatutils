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
| Weather | `process_weather_daymet`, `process_weather_gridmet`, `process_weather_nasapower`, `process_weather_openmeteo`, `process_weather_agera5`, `process_weather_nasapower_chirps`, `process_weather_dwd`, `process_weather_eobs`, `process_weather_xavier`, `process_weather_cmfd` |
| Soil | `process_soils_ssurgo`, `process_soils_gnatsgo`, `process_soils_isdasoil`, `process_soils_lucas`, `process_soils_soilgrids`, `process_soils_soilgrids_online`, `process_soils_hwsd` |

Coverage notes: Daymet = North America; GridMET/SSURGO/**gNATSGO** = USA;
DWD = Germany; E-OBS = Europe; **LUCAS** = EU; **Xavier (BR-DWGD)** = Brazil;
**CMFD** = China; **iSDAsoil** = Africa; NASA POWER / Open-Meteo / AgERA5 /
SoilGrids / HWSD2 = global. AgERA5 needs a free Copernicus CDS API key; CHIRPS
fuses NASA POWER with high-res rainfall (50S–50N).

**Regional sources added on top of the global set:**

| Source | Function | Region | Key / download |
|---|---|---|---|
| iSDAsoil | `process_soils_isdasoil` | Africa, 30 m | **none** — public S3 COGs streamed per-point |
| Xavier (BR-DWGD) | `process_weather_xavier` | Brazil, 0.1° | none, but **download** the NetCDFs (`xavier_nc_dir`) |
| CMFD | `process_weather_cmfd` | China, 0.1° 3-hourly | free **TPDC account** to download (`cmfd_nc_dir`) |
| LUCAS | `process_soils_lucas` | EU topsoil | free **ESDAC request form** to download (`lucas_csv`) |

`process_soils_isdasoil` mirrors the SSURGO signature (per-point `.SOL`, Saxton &
Rawls physics) using the iSDAsoil 0-20 / 20-50 cm predictions; `process_soils_lucas`
takes the nearest measured LUCAS topsoil sample (0-20 cm, **extrapolated** to the
rooting depth — flagged in the `.SOL`). `process_weather_xavier` reads the BR-DWGD
grids (solar radiation already in MJ/m²/day) and `process_weather_cmfd` aggregates
CMFD's 3-hourly fields to the daily statistics DSSAT needs (Tmax/Tmin, precip total,
mean radiation, RH/dew point from specific humidity).

**gNATSGO** (`process_soils_gnatsgo`) is the gap-free 30 m USDA grid (SSURGO +
STATSGO2 + Raster Soil Surveys): it returns a soil profile *everywhere there is
land*, filling the un-surveyed holes plain SSURGO leaves. The map-unit key is
read from the SoilWeb WCS grid; the tabular horizons come from the same USDA
Soil Data Access endpoint SSURGO uses, so where the underlying map unit matches,
the two produce byte-identical `.SOL` files. Same signature as
`process_soils_ssurgo`. No key required.

**DWD** (`process_weather_dwd`) downloads quality-controlled daily station
observations from the German Weather Service Open Data (no key) and uses the
nearest station with data for each point. Solar radiation is estimated from
sunshine duration via the Ångström–Prescott relation (FAO-56), since the daily
`kl` product does not measure it directly.

**E-OBS** (`process_weather_eobs`) reads the ECA&D 0.1° European gridded daily
NetCDFs. It includes daily global radiation (`qq`), so SRAD comes straight from
the data. Default mode points `eobs_nc_dir` at pre-downloaded E-OBS files (free
registration at www.ecad.eu); `eobs_use_cds=True` fetches an area subset via the
Copernicus CDS instead (needs the same `~/.cdsapirc` key as AgERA5).

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

### Optional backends (only for specific sources)

Most sources work with the base install. A few need an extra package:

| Source | R | Python |
|---|---|---|
| AgERA5, E-OBS (CDS mode) | `ecmwfr` (+ `~/.cdsapirc` key) | `dssatutils[agera5]` → `cdsapi` (+ key) |
| HWSD — SQLite DB | `DBI`, `RSQLite` | base (stdlib `sqlite3`) |
| HWSD — FAO Access `.mdb` directly | `odbc` (or `RODBC`) + the OS Microsoft Access ODBC driver | `dssatutils[hwsd]` → `pyodbc` + the Access ODBC driver |
| LUCAS `.xlsx` table | `readxl` | base (`pandas`/`openpyxl`) |

In a consumer repo that uses **renv**, these R packages are installed by the repo's
`setup_renv.R` (re-run it, or `renv::install(c("DBI","RSQLite","odbc","readxl")); renv::snapshot()`
to add them to an existing lockfile).

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
