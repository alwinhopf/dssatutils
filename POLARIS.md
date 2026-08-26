# POLARIS runtime and testing notes

The R POLARIS backend reads public POLARIS v1.0 GeoTIFF tiles through GDAL `/vsicurl` and `terra`. A sparse regional grid may touch many 1-degree tiles; each profile requires multiple property/depth rasters, so network responsiveness matters.

Runtime controls:

- `POLARIS_TIMEOUT_SEC` — GDAL HTTP timeout per remote read; default `45` seconds.
- `POLARIS_RETRIES` — application-level attempts per raster read; default `3`.
- `POLARIS_PROGRESS_EVERY` — emit progress after this many raster reads; default `10`.
- `cache_dir` / pipeline `polaris_cache_dir` — optionally persist downloaded GeoTIFFs rather than streaming them on every run.

The extractor groups points by POLARIS tile, reports the total number of raster reads before starting, reports each tile, and periodically reports completed reads. A raster that remains unavailable after all attempts contributes `NA` values instead of blocking the entire process indefinitely; profiles without usable texture are subsequently rejected.

Testing is split deliberately:

1. `tests/testthat/test_polaris.R` exercises tile addressing, transforms, van Genuchten/Saxton-Rawls water limits, `terra` extraction semantics using a local synthetic raster, bounded failure handling, runtime-setting validation, and DSSAT `.SOL` output.
2. `tests/test_polaris.py` covers the Python twin and checks key R/Python parity markers.
3. `tests/testthat/test_polaris_live.R` is an opt-in network integration test. With `RUN_POLARIS_LIVE=true`, it samples the public POLARIS median 0–5 cm clay raster at a known Iowa point through the same `/vsicurl` + `terra` path used in production.

The GitHub Actions smoke workflow runs the offline R suite on Linux, macOS, and Windows; the Python suite on the same three OS families with Python 3.9 and 3.12; and a separate Linux live-POLARIS R integration job.