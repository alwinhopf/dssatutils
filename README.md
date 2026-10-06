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
downloaded data. Both R and Python process one tile across all requested years,
write its complete point weather files, and release its assembled data before
starting the next tile. Download concurrency applies to years within that tile;
point-history memory therefore scales with points per tile and requested years,
rather than all points in the study. Existing annual CSV caches, completeness
checks, and atomic weather-file publication are preserved. A missing year blocks
publication for affected points while other tiles can still finish. Large tiles
or long histories can still need substantial memory.

The legacy `gridded` backend remains available for callers that need the original
daily-NetCDF ZIPs.

## Versioning

Semantic versioning with Git tags. **Consumer repos always pin to a tag**
(`@vX.Y.Z`), never `main`, so upstream changes never break a pipeline until you
deliberately bump the pin. Workflow: branch → CI smoke tests → merge → tag
`vX.Y.Z` → bump the pin in each consumer repo.

## Known limitations / notes

### Input-cache integrity fixes (2026-08-31)

- R GRIDMET extraction now retains an NA row for each out-of-coverage point.
  Previously `terra` dropped those cells and subsequent weather rows were
  assigned to the wrong point IDs. Python already retained invalid rows; its
  extraction helper now shares the same regression cases (invalid points at
  either end/in the middle, duplicate cells, out-of-order cells, single layers).
  Native NetCDF caches are reusable. Existing affected `.WTH` files are NOT
  repaired by installing this change: re-extract into a new weather directory,
  using the same native cache, then rebuild dependent model runs/results.
- SSURGO/gNATSGO's older R vector writer inserted an extra space before each
  layer after the first. Although the row writer was corrected on August 30,
  existing caches remained reusable under whitespace-only preflight checks.
  In DSSAT's fixed columns an intended 200 cm depth could become 20 cm.
- Alderman had a separate R/Python writer error: SLB occupied five rather than
  six columns. Full-width values could lose their first character; for example,
  SSKS 10.08 became 0.08. The writer now uses six-column SLB. All three USDA
  writers also use DSSAT's profile/site header widths; conductivity values of
  100 or more use one decimal to fit five numeric columns. Alderman missing
  conductivity/decimal values retain the `-99` sentinel without overflowing.
- `soil_file_issue(path)` returns `NULL`/`None` or a diagnostic, reading actual
  header-defined columns, not whitespace-separated apparent depths. It checks
  finite layer fields, alignment, positive/increasing depths and the 19-layer
  limit. This is a **format check**, not a physiological-quality guarantee.
  The gridded driver uses it before model execution. It does not change files.

**Offline recovery, without new downloads.** Keep original caches as evidence.
For SSURGO and gNATSGO, their existing layer mapping CSVs contain the properties
needed to regenerate `.SOL` files. Both languages expose the same helper:

```r
rebuild_soil_files_from_mapping("soil/study_SSURGO.CSV",
                              "soil/study_SSURGO_rebuilt", "SSURGO")
```

```python
from dssatutils import rebuild_soil_files_from_mapping
rebuild_soil_files_from_mapping("soil/study_SSURGO.CSV",
                              "soil/study_SSURGO_rebuilt", "SSURGO")
```

The destination must not exist (even as an empty directory). The result has
`ID` and `path` columns, preserves leading-zero IDs and stored property values,
and does not rerun pedotransfer calculations or download data. Use `GNATSGO`
for that source. **Alderman's mapping CSV contains metadata, not full layers**:
its legacy `.SOL` layer text must instead be read using the historical
five-column SLB layout and passed to the corrected Alderman writer. Do not
discard those files. Do not apply a generic formatter or whitespace parser to
them: horizon names can contain spaces.

R/Python parity note: validation, USDA writer corrections, recovery API and
test fixtures are mirrored. Real-cache verification rebuilt 143 SSURGO/gNATSGO
profiles identically byte-for-byte in R and Python; 72 Alderman profiles were
also reconstructed and validated in an isolated diagnostic. Native DSSAT
smoke runs completed with each of the three corrected soil writer outputs.
The 300/200 km sweep was not rerun and no production cache was overwritten.
Tests: `tests/test_input_cache_integrity.py` and its `tests/testthat/` twin.

Install/reload the corrected shared package in a **new** session after any
active sweep ends. Pinned releases do not yet contain uncommitted local fixes;
the current R session does not automatically reload edited package source.

- **Optional weather repair and QA** is available after provider downloads via
  `repair_weather_missing_values()`, `repair_weather_date_gaps()`,
  `repair_weather_temperature_inversions()`, and `audit_weather_quality()`.
  Repair functions support both neighbor averaging and bounded same-day swap
  (`method = "neighbor"` or `"swap"`, with `max_inversion_c = 2.0`);
  the audit writes flag-only findings to CSV and appends notes to
  `weather_repair.log`. Provider `NA`/`NaN`/infinite values are normalized to
  DSSAT's numeric `-99` marker before writing, and the repair reader also accepts
  those literal tokens in older cached files.
- **Weather-file validation** via `is_wth_valid()` understands DSSAT's
  fixed-width daily rows (including adjacent negative fields), requires
  consecutive dates, validates optional `start_year` and `end_year` bounds, and
  rejects physically impossible forcing while retaining the standard `-99` missing-value sentinel.
  Callers can pass `required_columns` to reject `-99` in model-essential fields while still
  permitting optional missing humidity/wind inputs.
- **AgERA5 time-series retrieval & assembly**:
  The time-series backend enforces per-cell CSV cache validation (`_valid_agera5_timeseries_csv`),
  atomic lock acquisition across workers, and inspection of the canonical CSV destination
  after requests return (even on non-path return objects). Incomplete annual jobs are not
  silently skipped; the assembler verifies exact calendar date coverage across the requested
  years before serialization. Incomplete points are reported without publishing partial
  `.WTH` files; the engine records exhausted attempts as retryable failures. Callers can enable `cache_only = True` for offline sweeps.
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

### Weather recovery integrity (2026-09-06)

R and Python now interpret `start_year` as January 1 and `end_year` as
December 31. `start_date` / `end_date` explicitly override those endpoints for
partial-year requests; a consecutive superset is allowed, but missing boundary
days are not. AgERA5 engine validation freezes an explicit endpoint capped at
its existing ten-day availability allowance. This does not change that allowance
or the configured temperature-repair policy.

Annual AgERA5 CSV checks reject any nonfinite or `-99` forcing value, invalid
coordinates, dates missing within a cell, duplicates, and cells outside the
requested area. Provider assembly also checks all seven forcing variables before
publication. Temperature inversions remain raw and are handled by optional shared
QC, not silently corrected by the downloader.

Time-series cache access uses OS byte-range locks compatible across R (`filelock`,
an R dependency) and Python (standard-library locking plus a thread guard).
A process exit releases ownership; elapsed time never steals a live lock. Empty
`.lock` files intentionally remain and must not be deleted during concurrent use.
Old directory-style `.lock` entries are not stolen: remove those only after
confirming their previous owners have stopped. Deploy/restart all workers together;
older package versions do not implement this protocol. Downloads use private
staging directories, validate before same-filesystem publication, and preserve
invalid prior cache bytes in `.invalid-*` evidence files when replacing them.
Cache-only misses do not delete prior files. R PSOCK workers explicitly resolve
exported downloader helpers rather than mixing old installed namespace functions.

Offline regressions: `test_weather_recovery_integrity.py` and its R twin cover
boundary/leap dates, individual missing values, multi-cell tiles, staged client
errors, cache preservation, and parallel assembly. Python tests additionally run
R/Python lock contention and process-exit recovery on the local platform.

### Offline CI compatibility

The offline suite supports Python 3.9 and newer; AgERA5 annotations are deferred
to avoid evaluating newer union syntax during Python 3.9 imports. This is a
Python compatibility detail with no R/Python behavioral divergence. R AgERA5 and
ERA5-Land conversion tests mock credential setup as well as downloads, so they
require neither CDS secrets nor a developer keyring.

## Daily NetCDF input contract

Generic gridded adapters require declared units and geographic, one-dimensional
latitude/longitude coordinates. Points outside the coordinate-centre domain,
projected grids, ambiguous variables, duplicate dates and subdaily records are
rejected. Reproject and aggregate upstream with the variable's physical semantics.
Radiation accepts J/kJ/MJ per square metre per daily record or daily-mean W/m²;
unknown units are errors. Vapour pressure accepts Pa/hPa/kPa/mbar. Wind must
declare its measurement height via the variable attribute `height_m` (metres)
or the adapter variable specification; m/s and km/h are supported and converted
to the WTH 2 m reference height. Already-MJ radiation is never rescaled.

SoilGrids online accepts `use_rest_api` explicitly in R and Python; an explicit
argument takes precedence over the legacy global/configuration default.

## Local and CI validation lanes

`scripts/pre-push.sh` runs the exact offline Python PR gate and records JUnit
and a test log. Native R tests run in the separate language-parity lane, with
logs uploaded on success or failure. Platform and live-provider checks remain
separate. Workflow actions are revision-pinned; Python 3.11 and R 4.3 are explicit.

## Wind units in DSSAT weather files

DSSAT v4.8 `WEATHER.CDE` specifies `WIND` in km/day. NASA POWER, its CHIRPS
hybrid, and AgERA5 provide wind in m/s; both language writers convert finite
nonnegative wind by 86.4 at serialization and preserve missing sentinels.
Weather validation and QA use km/day bounds (100 and 75 m/s equivalents).
Existing WTH caches from these adapters are not automatically rescaled: their
units must be verified against their writer provenance before regeneration or
a backed-up conversion, to avoid converting an already-correct file twice.
Raw provider caches remain in provider-native units. The October 6 correction
extends serialization-unit checks to all native weather writers with available wind.

Daymet R station metadata now uses the same column alignment as Python. An
extra leading space previously truncated coordinate precision in DSSAT
although whitespace-based validation saw the intended coordinate. Both
entry-point tests check the actual fixed-column coordinate fields. Previously
generated files require explicit header repair or regeneration.

The five active comparison adapters (Daymet, GridMET, NASA POWER, NASA POWER
+ CHIRPS v2, and AgERA5) reserve nine columns for station longitude in both
languages. This accommodates longitudes west of 100°W without shifting ELEV,
TAV or AMP across DSSAT's header boundaries. Regression fixtures include
three-digit negative longitudes; historical files require an explicit audit.
AgERA5 also retains two temperature decimals when both values would round to
zero at one decimal. This prevents rounding alone from creating a zero pair
that DSSAT rejects; true source zeros are preserved, not artificially perturbed.

GridMET AMP now uses the full mean annual range of monthly temperatures,
matching Daymet, POWER and AgERA5. DSSAT `SPAM/STEMP.for` divides TAMP by two
internally; the former `DSSAT::calc_AMP` convention halved it a second time.
Python and R share the corrected convention and an offline seasonal regression.
Previously cached headers require explicit regeneration.

### Fixed-column weather repair and macOS Alderman execution (2026-10-06)

R and Python weather QA, missing-value, date-gap and temperature-inversion
repairs now read the seven-character date and six-character numeric fields
before falling back to legacy whitespace rows. Adjacent values such as
`76.91062.7` represent RH2M 76.9 and WIND 1062.7 km/day when read.
Malformed rows raise an error instead of silently dropping days. Rewritten
large values now reserve the DSSAT separator and use integer precision;
parsing a legacy file alone does not fix the native model input.
All repair writers share the provider safeguard that retains two-decimal
near-freezing temperature pairs, so an unrelated repair cannot round them
back to a DSSAT-rejected pair of zeros. Trailing whitespace remains supported.
The October 5 wind-unit correction exposed whitespace reader assumptions dating
to June; the wind conversion is retained.

On macOS, the R SSURGO, gNATSGO and Alderman adapters execute serially even
when multiple cores are requested, avoiding the fork path after native libraries initialize.
Windows retains socket workers and Linux retains its existing parallel path.
Python uses its existing execution path; this is an R-specific platform
implementation difference, not a soil-output schema change. Refresh the local
R package and restart R before rerunning failed scenarios; editing source does
not update an already loaded package.

## DSSAT serialization and platform corrections (6 October 2026)

Daily numeric fields preserve one separator and a five-character token, matching
DSSAT's header-based column reader. Writers prefer one decimal (two for the
existing near-freezing safeguard), reducing precision when rounding would consume
the separator: 999.96 becomes 1000 and 1062.7 becomes 1063. Negative temperatures
and dewpoints remain valid. Nonfinite values use -99; unrepresentable values raise
instead of clipping or shifting columns. The neutral writer also retains the
four-digit-year `$WEATHER` marker and seven-character date header.

Provider data and unit/height conversion helpers retain m/s internally. All
native DSSAT writers with available wind convert finite nonnegative m/s to
km/day exactly once (multiply by 86.4); missing wind remains -99. This includes
Open-Meteo, ERA5-Land, CMFD, DWD, E-OBS, Xavier, the common gridded writer and
NASA POWER rainfall hybrids. Height adjustments retain matching WNDHT metadata;
ERA5-Land explicitly declares its 10 m wind measurement height. Existing provider
caches are not rescaled automatically: verify their writer provenance before
regenerating outputs, so correctly converted wind is never converted twice.

QA and repairs recover adjacent legacy six-character fields without dropping
days. Rewritten files reserve the separator expected by the native reader;
large values consequently have integer precision. Existing adjacent-field WTH
files require reserialization before native model use; parsing alone does not
correct the file consumed by DSSAT. iSDAsoil, LUCAS, SSURGO, gNATSGO and Alderman
conductivity writers also handle the rounding boundary at 100 and character-valued
numeric responses. Existing provider conductivity caps are unchanged.

The R SSURGO/gNATSGO/Alderman adapters use sequential work on macOS, socket
workers on Windows and the existing fork path on Linux. This R-specific platform
policy preserves the R/Python output contract. Offline regressions exercise
native header columns, legacy row recovery, missing markers, conductivity and
socket workers with shared fixtures in both languages.

Native R parity tests now run separately on Linux, Windows and macOS, alongside
the Python platform jobs. CI passes the interpreter's native absolute path from
Python into R (including Windows) and uploads R JUnit XML and logs on failures.
Pytest package gates prioritize their own checkout's source to avoid testing a
different editable installation. Hosted results must be checked after publishing
changes; local macOS validation does not establish hosted Windows/Linux success.

The combined October 6 results and hosted-CI limits are recorded in
[the monorepo validation report](../cropmodel/WEATHER_SOIL_PORTABILITY_VALIDATION.md).
