# AGENTS.md — dssatutils

> **Workspace context:** Read the root [`../AGENTS.md`](../AGENTS.md) first. This document
> holds rules and guidance specific to the `dssatutils` repository.

## 1. Role in the Workspace

`dssatutils` is the **foundation weather & soil download/conversion library** for DSSAT:
- It downloads gridded, station, and reanalysis data from 25+ weather sources and 15+ soil databases.
- It formats extracted profiles into valid DSSAT `.WTH` (weather) and `.SOL` (soil) files.
- It provides weather quality auditing, date gap filling, and missing value repair.
- It is consumed by `DSSAT_Gridded_Run_Tutorial`, `dssatcalibrator`, `DSSAT_ML_Phenology_Prediction`, and `DSSAT_SubField_MILP_Analysis`.

## 2. 1:1 R ↔ Python Parity Contract

`dssatutils` is a dual-language package where every public source adapter exists in both languages:
- R implementation: `R/weather_*.R`, `R/soil_*.R`, `R/credentials.R`, `R/config.R`.
- Python implementation: `python/dssatutils/weather_*.py`, `python/dssatutils/soil_*.py`, `python/dssatutils/credentials.py`, `python/dssatutils/config.py`.

### When adding or editing a data source, update all 6 components together:
1. R implementation under `dssatutils/R/`.
2. Python implementation under `dssatutils/python/dssatutils/`.
3. `dssatutils/NAMESPACE` (`export(process_...)`).
4. `dssatutils/python/dssatutils/__init__.py` (`_EXPORTS` map and `__all__`).
5. Offline tests in `tests/testthat/` (R) and `tests/` (Python).
6. Documentation in `README.md` and tutorial `DATA_SOURCES.md`.

## 3. Critical Implementation Rules & Common Pitfalls

### A. Fixed-Width DSSAT File Formatting
- DSSAT `.WTH` and `.SOL` files use strict fixed-width Fortran column positions.
- **Never reformat `.WTH` or `.SOL` files using generic code formatters or lashing whitespace tools.**
- Soil layers must have matching depths, layer counts, and valid physical ranges (e.g. `SLDP`, `SLLL`, `SDUL`, `SSAT`).

### B. Python PEP 562 Lazy Loading
- Python submodules pull in heavy GIS packages (xarray, rasterio, geopandas, netCDF4, cdsapi) only when imported.
- `python/dssatutils/__init__.py` uses `__getattr__` to load submodules on-demand. Never eagerly import heavy GIS modules at package root level.

### C. Offline Testing with Committed Fixtures
- All automated tests must run offline without live network access or real API keys.
- Mock API responses and NetCDF slices using committed fixtures in `tests/fixtures/`.
- Live download tests must be marked with `@pytest.mark.live` or gated behind environment flags.

### D. Credentials & CDS Setup
- Never commit API keys, tokens, or personal paths.
- CDS-backed sources (AgERA5, ERA5-Land, E-OBS CDS mode) configure authentication via `setup_cds_credentials()`, reading `CDSAPI_KEY`/`CDSAPI_URL` or `~/.cdsapirc`.

## 4. Verification & Testing

```bash
# Python tests
python -m pytest tests/test_new_sources.py tests/test_smoke.py -q
python -m pytest tests/test_comprehensive.py -q

# R tests
Rscript -e "pkgload::load_all('.', quiet=TRUE); testthat::test_file('tests/testthat/test_new_sources.R')"
Rscript -e "for (f in list.files('R', pattern='\\.R$', full.names=TRUE)) parse(f)"
```
